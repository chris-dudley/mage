#include "yens.hpp"

#include <algorithm>
#include <queue>
#include <functional>

namespace shortest_paths {

bool greater_cost( const Path<>& lhs, const Path<>& rhs ) {
    return lhs.total_cost > rhs.total_cost;
}

// Check if our path heap contains the given path.
bool contains_path( const std::vector<Path<>>& paths, const Path<>& to_find ) {
    for (const auto& path : paths) {
        if (path.edges == to_find.edges) {
            return true;
        }
    } 

    return false;
}

// Returns the indicies of paths in `paths` that have `total_cost == target_cost`, up to `limit` entries.
std::vector<size_t> find_same_cost_paths(const std::vector<Path<>>& paths, double target_cost, uint64_t limit) {
    std::vector<size_t> result;
    for (size_t i = 0; i < paths.size() && result.size() < limit; i++) {
        if (paths[i].total_cost == target_cost) {
            result.push_back(i);
        }
    }
    return result;
}

std::vector<Path<>> KShortestPaths(
    const mg_graph::GraphView<> &graph, std::uint64_t source_id, std::uint64_t sink_id,
    std::uint64_t K, ShortestPathFunc shortest_path_func,
    CheckAbortFunc check_abort
) {
    std::vector<Path<>> result;
    // Using a heap instead of a priority_queue so that we can search the heap for
    // duplicate paths.
    std::vector<Path<>> possible_paths;

    // Return early in cases with no paths
    if (source_id == sink_id || K == 0) {
        return result;
    }

    auto shortest_path = shortest_path_func(graph, source_id, sink_id, {}, {}, check_abort);
    if (shortest_path.empty()) {
        // No paths found
        return result;
    }
    result.push_back(shortest_path);

    if (K == 1) {
        // Don't bother with the rest of the algorithm if K == 1
        return result;
    }
    
    for (uint64_t k = 1; k < K; k++) {
        // k-1'th path
        const auto& prev_shortest = result.back();
        uint64_t paths_still_needed = K - k;

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
            EdgeIdSet ignored_edges;
            NodeIdSet ignored_nodes;

            // Find all previous shortest paths that share same root path edges and remove the edge
            // used to go to the next node in that path.
            for (const auto& prev_path : result) {
                if (prev_path.has_prefix(root_path)) {
                    // Ignore outgoing edge from spur node to next node in previous shortest path
                    ignored_edges.insert(prev_path.edges[spur_index]);
                }
            }

            // Ignore all nodes in the root path except the spur node
            for (size_t i = 0; i < spur_index; i++) {
                ignored_nodes.insert(root_path.nodes[i]);
            }

            // Try all remaining spurs from this node
            size_t remaining_out_edges = 0;
            for (const auto& neighbor : graph.OutNeighbours(spur_node)) {
                if (!ignored_edges.contains(neighbor.edge_id)) {
                    remaining_out_edges += 1;
                }
            }
            if (remaining_out_edges == 0) {
                // All edges out from this node exhausted, continue
                continue;
            } 

            while (remaining_out_edges-- > 0) {
                auto spur_path = shortest_path_func(graph, spur_node, sink_id, ignored_edges, ignored_nodes, check_abort);
                if (spur_path.empty()) {
                    // No more available paths from this node.
                    break;
                }

                // Ignore newly taken edge for further attempts
                ignored_edges.insert(spur_path.edges[0]);

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

} // namespace shortest_paths