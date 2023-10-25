#pragma once

#include <functional>
#include <algorithm>
#include <queue>
#include <thread>
#include <memory>

#include <omp.h>

#include "shortest_path.hpp"
#include "dijkstra.hpp"

namespace shortest_paths {

/// @brief Class implementing Yen's K-shortest paths algorithm.
/// @tparam TSize Type used for node and edge IDs.
/// @tparam Pathfinder Pathfinder used to find shortest path between a given source and target.
template<typename TSize = std::uint64_t, typename Pathfinder = DijkstraPathfinder<TSize>>
class YensPathfinder {
public:
    /// @brief Type of graph view this pathfinder expects.
    using GraphViewType = mg_graph::GraphView<TSize>;
    using Self = DijkstraPathfinder<TSize>;
    /// @brief Represents a set of edges by ID.
    using EdgeIdSet = std::unordered_set<TSize>;
    /// @brief Represents a set of nodes by ID.
    using NodeIdSet = std::unordered_set<TSize>;
    /// @brief A vector of nodes by ID.
    using NodeIdVec = std::vector<TSize>;
    /// @brief A vector of edges by ID.
    using EdgeIdVec = std::vector<TSize>;
    /// @brief A vector of paths.
    using PathVec = std::vector<Path<TSize>>;

    /// @brief Search for the K-shortest paths.
    /// @param graph The Graph to search.
    /// @param source_id ID of source node for pathfinding.
    /// @param target_id ID of target node for pathfinding.
    /// @param K Number of shortest paths to search for.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @param threads Maximum number of threads to use during pathfinding.
    /// @return A vector containing the shortest paths in order of increasing total cost.
    PathVec search(
        const GraphViewType& graph, TSize source_id, TSize target_id, size_t K,
        const CheckAbortFunc& check_abort = CheckAbortNoop, int threads = 0
    ) {
        EdgeIdSet ignored_edges;
        NodeIdSet ignored_nodes;
        return search(graph, source_id, target_id, K, ignored_edges, ignored_nodes, check_abort, threads);
    }
    
    /// @brief Search for the K-shortest paths.
    /// @param graph The Graph to search.
    /// @param source_id ID of source node for pathfinding.
    /// @param target_id ID of target node for pathfinding.
    /// @param K Number of shortest paths to search for.
    /// @param ignored_edges IDs of edges to ignore during pathfinding.
    /// @param ignored_nodes IDs of nodes to ignore during pathfinding.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @param threads Maximum number of threads to use during pathfinding.
    /// @return A vector containing the shortest paths in order of increasing total cost.
    PathVec search(
        const GraphViewType& graph, TSize source_id, TSize target_id, size_t K,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort = CheckAbortNoop, int threads = 0
    ) {
        if (threads == 1) {
            return search_sequential(graph, source_id, target_id, K, ignored_edges, ignored_nodes, check_abort);
        } else {
            return search_omp(graph, source_id, target_id, K, ignored_edges, ignored_nodes, check_abort, threads);
        }
    }

private:
    using PathfinderUniquePtr = std::unique_ptr<Pathfinder>;

    struct Task {
        TSize spur_id;
        Path<TSize> root_path;
        EdgeIdSet ignored_edges;
        NodeIdSet ignored_nodes;
    };

    PathVec search_sequential(
        const GraphViewType& graph, TSize source_id, TSize target_id, size_t K,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort
    ) {
        PathVec result;
        // Using a heap instead of a priority_queue so that we can search the heap for
        // duplicate paths.
        PathVec possible_paths;

        // Return early in cases with no paths
        if (source_id == target_id || K == 0) {
            return result;
        }

        Pathfinder pathfinder;
        pathfinder.search(graph, source_id, target_id, ignored_edges, ignored_nodes, check_abort);
        if (!pathfinder.has_path_to(target_id)) {
            // No paths founds
            return result;
        }
        auto shortest_path = pathfinder.path_to(target_id);
        result.push_back(shortest_path);

        if (K == 1) {
            // Don't bother with the rest of the algorithm if K == 1
            return result;
        }

        for (uint64_t k = 1; k < K; k++) {
            // k-1'th path
            const auto& prev_shortest = result.back();
            size_t paths_still_needed = K - k;

            // Check if we have at least `paths_still_needed` paths with the same weight as `prev_shortest`
            // in `possible_paths`. If so, we will not find any shorter paths than that, so just add those
            // to the result and return.
            if (possible_paths.size() >= paths_still_needed) {
                auto same_cost_path_indicies = find_same_cost_paths(possible_paths, prev_shortest.total_cost, paths_still_needed);
                if (same_cost_path_indicies.size() == paths_still_needed) {
                    for (size_t i = 0; i < same_cost_path_indicies.size(); i++) {
                        result.push_back(std::move(possible_paths[i]));
                    }
                    return result;
                }
            }

            // The spur node ranges from the first node to the next to last node in the previous
            // k-shortest path.
            
            // We can skip any edges that are shared with the k-2'th path, since we've already
            // calculated spurs from those source nodes. (Lawler's modification).
            uint64_t spur_index = 0;
            if (k > 1) {
                const auto& ancestor = result[k-2];
                const auto max_prefix_size = std::min(ancestor.size(), prev_shortest.size());
                while (
                    spur_index < max_prefix_size
                    && prev_shortest.edges[spur_index] == ancestor.edges[spur_index]
                ) {
                    spur_index++;
                }
            }
            for (; spur_index < prev_shortest.size(); spur_index++) {
                // Check if we should abort before each shortest path call, as it could take a while
                check_abort();

                // Note: On this path, we are not taking the edge at prev_shortest[spur_index]
                const auto spur_node = prev_shortest.nodes[spur_index];
                // Re-use spur_index for the number of edges we want to keep.
                // If spur is index 0, we want 0 previous edges, if index 1 we want 1 edge, etc.
                auto root_path = prev_shortest.prefix(spur_index);

                // Ignoring edges and nodes instead of removing them from the graph to avoid having to
                // rebuild the entire graph view every time, as you can't copy the view.
                EdgeIdSet cur_ignored_edges(ignored_edges);
                NodeIdSet cur_ignored_nodes(ignored_nodes);

                // Find all previous shortest paths that share same root path edges and remove the edge
                // used to go to the next node in that path.
                for (const auto& prev_path : result) {
                    if (prev_path.has_prefix(root_path)) {
                        // Ignore outgoing edge from spur node to next node in previous shortest path
                        cur_ignored_edges.insert(prev_path.edges[spur_index]);
                    }
                }

                // Ignore all nodes in the root path except the spur node
                for (size_t i = 0; i < spur_index; i++) {
                    cur_ignored_nodes.insert(root_path.nodes[i]);
                }

                // Try all remaining spurs from this node
                size_t remaining_out_edges = 0;
                for (const auto& neighbor : graph.OutNeighbours(spur_node)) {
                    if (!cur_ignored_edges.contains(neighbor.edge_id)) {
                        remaining_out_edges += 1;
                    }
                }
                if (remaining_out_edges == 0) {
                    // All edges out from this node exhausted, continue
                    continue;
                } 

                while (remaining_out_edges-- > 0) {
                    pathfinder.search(graph, spur_node, target_id, cur_ignored_edges, cur_ignored_nodes, check_abort);
                    auto spur_path = pathfinder.path_to(target_id);
                    if (spur_path.empty()) {
                        // No more available paths from this node.
                        break;
                    }

                    // Ignore newly taken edge for further attempts
                    cur_ignored_edges.insert(spur_path.edges[0]);

                    // Join together the root and spur paths
                    auto total_path = root_path.join(spur_path);

                    // If total_path is not in possible_paths, add it to the heap
                    if (!contains_path(possible_paths, total_path)) {
                        possible_paths.push_back(total_path);
                        std::push_heap(possible_paths.begin(), possible_paths.end(), greater_cost);
                    }
                }
            }

            if (possible_paths.empty()) {
                // No spur paths left, don't bother checking for any more paths.
                break;
            }

            // Pop the minimum weighted path from the heap and add it to the result
            std::pop_heap(possible_paths.begin(), possible_paths.end(), greater_cost);
            result.push_back(std::move(possible_paths.back()));
            possible_paths.pop_back();
        }

        return result;
    }

    PathVec search_omp(
        const GraphViewType& graph, TSize source_id, TSize target_id, size_t K,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort, int threads
    ) {
        if (threads <= 0) {
            threads = std::thread::hardware_concurrency();
        }

        PathVec result;
        // Using a heap instead of a priority_queue so that we can search the heap for
        // duplicate paths.
        PathVec possible_paths;

        // Return early in cases with no paths
        if (source_id == target_id || K == 0) {
            return result;
        }

        // Set up a vector of unique_ptr to pathfinders so we can re-use them later and reduce the amount
        // of memory allocation/freeing needed.
        std::vector<PathfinderUniquePtr> pathfinders;
        pathfinders.emplace_back(std::move(std::make_unique<Pathfinder>()));

        pathfinders[0]->search(graph, source_id, target_id, ignored_edges, ignored_nodes, check_abort);
        if (!pathfinders[0]->has_path_to(target_id)) {
            // No paths founds
            return result;
        }
        auto shortest_path = pathfinders[0]->path_to(target_id);
        result.push_back(shortest_path);

        if (K == 1) {
            // Don't bother with the rest of the algorithm if K == 1
            return result;
        }

        omp_set_dynamic(threads);

        for (uint64_t k = 1; k < K; k++) {
            // k-1'th path
            const auto& prev_shortest = result.back();
            size_t paths_still_needed = K - k;

            // Check if we have at least `paths_still_needed` paths with the same weight as `prev_shortest`
            // in `possible_paths`. If so, we will not find any shorter paths than that, so just add those
            // to the result and return.
            if (possible_paths.size() >= paths_still_needed) {
                auto same_cost_path_indicies = find_same_cost_paths(possible_paths, prev_shortest.total_cost, paths_still_needed);
                if (same_cost_path_indicies.size() == paths_still_needed) {
                    for (size_t i = 0; i < same_cost_path_indicies.size(); i++) {
                        result.push_back(std::move(possible_paths[i]));
                    }
                    return result;
                }
            }

            // Put together all the tasks for this round, as they only depend on the previous shortest paths
            // computed for previous values of k.
            std::vector<Task> tasks;

            // The spur node ranges from the first node to the next to last node in the previous
            // k-shortest path.
            
            // We can skip any edges that are shared with the k-2'th path, since we've already
            // calculated spurs from those source nodes. (Lawler's modification).
            uint64_t spur_index = 0;
            if (k > 1) {
                const auto& ancestor = result[k-2];
                const auto max_prefix_size = std::min(ancestor.size(), prev_shortest.size());
                while (
                    spur_index < max_prefix_size
                    && prev_shortest.edges[spur_index] == ancestor.edges[spur_index]
                ) {
                    spur_index++;
                }
            }
            for (; spur_index < prev_shortest.size(); spur_index++) {
                // Check if we should abort before each shortest path call, as it could take a while
                check_abort();

                // Note: On this path, we are not taking the edge at prev_shortest[spur_index]
                const auto spur_node = prev_shortest.nodes[spur_index];
                // Re-use spur_index for the number of edges we want to keep.
                // If spur is index 0, we want 0 previous edges, if index 1 we want 1 edge, etc.
                auto root_path = prev_shortest.prefix(spur_index);

                // Ignoring edges and nodes instead of removing them from the graph to avoid having to
                // rebuild the entire graph view every time, as you can't copy the view.
                EdgeIdSet cur_ignored_edges(ignored_edges);
                NodeIdSet cur_ignored_nodes(ignored_nodes);

                // Find all previous shortest paths that share same root path edges and remove the edge
                // used to go to the next node in that path.
                for (const auto& prev_path : result) {
                    if (prev_path.has_prefix(root_path)) {
                        // Ignore outgoing edge from spur node to next node in previous shortest path
                        cur_ignored_edges.insert(prev_path.edges[spur_index]);
                    }
                }

                // Ignore all nodes in the root path except the spur node
                for (size_t i = 0; i < spur_index; i++) {
                    cur_ignored_nodes.insert(root_path.nodes[i]);
                }

                // Try all remaining spurs from this node
                for (const auto& neighbor : graph.OutNeighbours(spur_node)) {
                    if (!cur_ignored_edges.contains(neighbor.edge_id)) {
                        // Create a copy of ignored_edges with all other outbound edges from this
                        // node removed so we try this specific edge.
                        EdgeIdSet ignored_edges_copy(cur_ignored_edges);
                        for (const auto& other_neighbor : graph.OutNeighbours(spur_node)) {
                            if (other_neighbor.edge_id != neighbor.edge_id) {
                                ignored_edges_copy.emplace(other_neighbor.edge_id);
                            }
                        }

                        tasks.emplace_back(spur_node, root_path, std::move(ignored_edges_copy), cur_ignored_nodes);
                    }
                }
            }

            // Add more pathfinders if needed
            for (size_t i = pathfinders.size(); i < tasks.size(); i++) {
                pathfinders.emplace_back(std::move(std::make_unique<Pathfinder>()));
            }

            // Perform pathfinding for all tasks
            if (tasks.size() > 1) {
                #pragma omp parallel for
                for (size_t task_id = 0; task_id < tasks.size(); task_id++) {
                    const auto& task = tasks[task_id];
                    pathfinders[task_id]->search(graph, task.spur_id, target_id, task.ignored_edges, task.ignored_nodes, check_abort);
                }
            } else if (tasks.size() == 1) {
                // Don't bother with overhead of an omp parallel for if only 1 task
                const auto& task = tasks[0];
                pathfinders[0]->search(graph, task.spur_id, target_id, task.ignored_edges, task.ignored_nodes, check_abort);
            }

            // Process results for all tasks
            for (size_t task_id = 0; task_id < tasks.size(); task_id++) {
                if (!(pathfinders[task_id]->has_path_to(target_id))) {
                    continue;
                }
                auto spur_path = pathfinders[task_id]->path_to(target_id);
                auto total_path = tasks[task_id].root_path.join(spur_path);

                if (!contains_path(possible_paths, total_path)) {
                    possible_paths.push_back(total_path);
                    std::push_heap(possible_paths.begin(), possible_paths.end(), greater_cost);
                }
            }

            if (possible_paths.empty()) {
                // No spur paths left, don't bother checking for any more paths.
                break;
            }

            // Pop the minimum weighted path from the heap and add it to the result
            std::pop_heap(possible_paths.begin(), possible_paths.end(), greater_cost);
            result.push_back(std::move(possible_paths.back()));
            possible_paths.pop_back();
        }

        return result;
    }

    // Returns true if lhs has a greater total cost than rhs
    static bool greater_cost(const Path<TSize>& lhs, const Path<TSize>& rhs) {
        return lhs.total_cost > rhs.total_cost;
    }

    // Returns the indicies of paths in `paths` that have `total_cost == target_cost`, up to `limit` entries.
    static std::vector<size_t> find_same_cost_paths(const PathVec& paths, double target_cost, size_t limit) {
        std::vector<size_t> result;
        for (size_t i = 0; i < paths.size() && result.size() < limit; i++) {
            if (paths[i].total_cost == target_cost) {
                result.push_back(i);
            }
        }
        return result;
    }

    // Check if our path heap contains the given path.
    static bool contains_path( const PathVec& paths, const Path<TSize>& to_find ) {
        for (const auto& path : paths) {
            if (path.edges == to_find.edges) {
                return true;
            }
        } 

        return false;
    };
};

/// @brief Computes the K shortest paths in the graph from source to sink.
/// @param graph Current graph.
/// @param source_id ID of source node for paths.
/// @param sink_id ID of final node for paths.
/// @param K Number of shortest paths to compute.
/// @param check_abort Function used to check if execution should be aborted.
/// @param threads Number of threads to use during pathfinding.
/// @return A vector of paths, each contaning edges taken in the path from source to sink.
std::vector<Path<>> KShortestPaths(
    const mg_graph::GraphView<> &graph, std::uint64_t source_id, std::uint64_t sink_id,
    std::uint64_t K, const CheckAbortFunc& check_abort, int threads = 0
);

} // namespace shortest_paths