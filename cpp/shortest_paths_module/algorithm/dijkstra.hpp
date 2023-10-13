#pragma once

#include <optional>
#include <tuple>
#include <limits>
#include <ranges>
#include <stdexcept>

#include <mg_graph.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include "shortest_path.hpp"

namespace shortest_paths {

/// @brief Class implementing Dijkstra's algorithm for computing single-source shortest paths.
/// @tparam TSize Type used for node and edge IDs.
template<typename TSize = std::uint64_t>
class DijkstraPathfinder {
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

    /// @brief Infinite distance used for unreachable nodes.
    static constexpr const double POSITIVE_INFINITY = std::numeric_limits<double>::infinity();
    /// @brief Amount distances must differ by in order to be considered different. Used to handle
    /// numerical instability that is inherent to floating point numerical representations.
    static constexpr const double EPSILON = 1e-14;
    /// @brief Marker for an invalid node in the parents vector.
    static constexpr const TSize INVALID_NODE = std::numeric_limits<TSize>::max();
    /// @brief Marker for an invalid edge in the vector of edges leading into nodes in the spanning tree.
    static constexpr const TSize INVALID_EDGE = std::numeric_limits<TSize>::max();

private:
    struct HeapData {
        TSize node_id;
        double distance;

        std::partial_ordering operator<=>(const HeapData& rhs) const {
            return distance <=> rhs.distance;
        }
    };

    // Type used for queue of nodes to visit, sorted by distance from source.
    // Making the comparison function std::greater makes it a min-queue.
    using HeapType = boost::heap::fibonacci_heap<HeapData, boost::heap::compare<std::greater<HeapData>>>;

    // `dist_to[i]` holds the shortest distance from `source` to `i`.
    std::vector<double> dist_to;
    // `added[i]` will be true if vertex `i` is included in the shortest path tree or
    // shortest distance from `source` to `i` is finalized.
    std::vector<bool> added;
    // `parent[i]` is the ID of the parent node of vertex `i` in the shortest path tree.
    std::vector<TSize> parent;
    // `edge_into[i]` is the ID of the edge used to reach node `i` in the shortest path tree.
    std::vector<TSize> edge_into;
    // ID of the source node for last pathfinding.
    TSize source_id;
    // ID of target node, if any.
    std::optional<TSize> target_id;

public:
    DijkstraPathfinder() = default;
    DijkstraPathfinder(const Self&) = default;
    DijkstraPathfinder(Self&&) = default;
    Self& operator=(const Self&) = default;
    Self& operator=(Self&&) = default;

    /// @brief Resets the state of the pathfinder and then searches for the shorest paths on `graph` from `source`
    ///     to all other reachable nodes.
    /// @param graph The graph to search.
    /// @param source The ID of the source node.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @throws std::invalid_argument if the source node is not in the graph.
    void search(const GraphViewType& graph, TSize source, const CheckAbortFunc &check_abort = CheckAbortNoop) {
        EdgeIdSet empty_edges;
        NodeIdSet empty_nodes;
        search(graph, source, empty_edges, empty_nodes, check_abort);
    }

    /// @brief Resets the state of the pathfinder and then searches for the shorest paths on `graph` from `source`
    ///     to all other reachable nodes.
    /// @param graph The graph to search.
    /// @param source The ID of the source node.
    /// @param ignored_edges Set of edge IDs to ignore during pathfinding.
    /// @param ignored_nodes Set of node IDs to ignore during pathfinding.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @throws std::invalid_argument if the source node is not in the graph.
    void search(
        const GraphViewType& graph, TSize source,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc &check_abort = CheckAbortNoop
    ) {
        TSize num_nodes = graph.Nodes().size();
        reset(num_nodes, source);

        if (num_nodes == 0) {
            // Empty graph, no paths to find
            return;
        }

        if (source >= num_nodes) {
            throw std::invalid_argument("source node not in graph");
        }
        do_search(graph, ignored_edges, ignored_nodes, check_abort);
    }

    /// @brief Resets the state of the pathfinder and then searches for the shorest paths on `graph` from `source`
    ///     to `target`, if such a path exists.
    /// @param graph The graph to search.
    /// @param source The ID of the source node.
    /// @param target The ID of the target node.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @throws std::invalid_argument if the source or target nodes are not in the graph.
    void search(const GraphViewType& graph, TSize source, TSize target, const CheckAbortFunc &check_abort = CheckAbortNoop) {
        EdgeIdSet empty_edges;
        NodeIdSet empty_nodes;
        search(graph, source, target, empty_edges, empty_nodes);
    }

    /// @brief Resets the state of the pathfinder and then searches for the shorest paths on `graph` from `source`
    ///     to `target`, if such a path exists.
    /// @param graph The graph to search.
    /// @param source The ID of the source node.
    /// @param target The ID of the target node.
    /// @param ignored_edges Set of edge IDs to ignore during pathfinding.
    /// @param ignored_nodes Set of node IDs to ignore during pathfinding.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @throws std::invalid_argument if the source or target nodes are not in the graph.
    void search(
        const GraphViewType& graph, TSize source, TSize target,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc &check_abort = CheckAbortNoop
    ) {
        TSize num_nodes = graph.Nodes().size();
        reset(num_nodes, source, target);

        if (num_nodes == 0) {
            // Empty graph, no paths to find
            return;
        }

        if (source >= num_nodes) {
            throw std::invalid_argument("source node not in graph");
        }
        if (target >= num_nodes) {
            throw std::invalid_argument("target not node in graph");
        }
        do_search(graph, ignored_edges, ignored_nodes, check_abort);
    }

    /// @brief Returns the number of reachable nodes in the graph, including the source.
    size_t num_reachable_nodes() const {
        size_t result = 0;
        for (size_t node_id = 0; node_id < dist_to.size(); node_id++) {
            if (dist_to[node_id] < POSITIVE_INFINITY) {
                result++;
            }
        }
        return result;
    }

    /// @brief Returns the set of nodes that are the source or reachable from the source.
    NodeIdVec reachable_nodes() const {
        std::vector<TSize> result;
        for (size_t node_id = 0; node_id < dist_to.size(); node_id++) {
            if (dist_to[node_id] < POSITIVE_INFINITY) {
                result.push_back(node_id);
            }
        }
        return result;
    }

    /// @brief Returns the set of edges used to traverse the minimum-weighted paths to reachable nodes.
    /// If the graph has a negative cycle, the returned edges are not guaranteed to contain all edges needed
    /// to reach all nodes.
    EdgeIdVec edges_used() const {
        EdgeIdVec result;
        for (const auto& maybe_edge : edge_into) {
            if (maybe_edge != INVALID_EDGE) {
                result.push_back(maybe_edge);
            }
        }
        return result;
    }

    /// @brief Returns a vector containing the distances from the source to other nodes in the graph.
    /// If a target was specified, not all reachable nodes may be populated.
    const std::vector<double>& distances() const {
        return dist_to;
    }
    /// @brief Returns whether a path from the source to the specified node exists.
    /// @param node The ID of the destination node.
    bool has_path_to(TSize node) const {
        if (node >= dist_to.size()) {
            return false;
        }
        return dist_to[node] < POSITIVE_INFINITY;
    }

    /// @brief Returns the lowest-cost path from the source to the specified node, if one exists.
    /// @param node The destination node.
    Path<TSize> path_to(TSize node) const {
        Path<TSize> result(source_id);
        if (!has_path_to(node) || node == source_id) {
            return result;
        }

        // (edge id, from, to, distance)
        using EdgeInfo = std::tuple<TSize, TSize, TSize, double>;
        std::vector<EdgeInfo> stack;

        // Build a stack of edges starting at the destination, then reverse it to build the path.
        for (TSize current_node = node; current_node != source_id; current_node = parent[current_node]) {
            stack.emplace_back(edge_into[current_node], parent[current_node], current_node, dist_to[current_node]);
        }

        for (auto [edge_id, from_node, to_node, distance] : std::views::reverse(stack)) {
            result.add_edge(edge_id, from_node, to_node, distance);
        }

        return result;
    }

private:
    void reset(TSize num_nodes, TSize source, std::optional<TSize> target = std::nullopt) {
        dist_to.assign(num_nodes, POSITIVE_INFINITY);
        added.assign(num_nodes, false);
        parent.assign(num_nodes, INVALID_NODE);
        edge_into.assign(num_nodes, INVALID_EDGE);
        source_id = source;
        target_id = target;
    }

    void do_search(
        const GraphViewType& graph,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc &check_abort
    ) {
        dist_to[source_id] = 0.0;

        HeapType node_queue;
        // Seed queue with source node
        node_queue.emplace(source_id, 0.0);

        while (!node_queue.empty()) {
            check_abort();

            auto current_node = node_queue.top().node_id;
            node_queue.pop();
            if (added[current_node]) {
                // Don't revisit nodes
                continue;
            }
            added[current_node] = true;

            if (target_id.has_value() && target_id.value() == current_node) {
                // Found the target, can stop early
                break;
            }

            relax(graph, current_node, node_queue, ignored_edges, ignored_nodes);
        }
    }

    void relax(
        const GraphViewType& graph, TSize current_node, HeapType& node_queue,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes
    ) {
        for (const auto& neighbor : graph.OutNeighbours(current_node)) {
            // Don't update parent + edge in if we've already added the neighbor to the shortest path,
            // or we should ignore the edge or neighbor node.
            if (added[neighbor.node_id] || ignored_edges.contains(neighbor.edge_id) || ignored_nodes.contains(neighbor.node_id)) {
                continue;
            }
            auto next_node = neighbor.node_id;
            double edge_weight = graph.IsWeighted() ? graph.GetWeight(neighbor.edge_id) : 1.0;
            double dist_to_next = dist_to[current_node] + edge_weight;

            if (dist_to[next_node] > (dist_to_next + EPSILON)) {
                dist_to[next_node] = dist_to_next;

                parent[next_node] = current_node;
                edge_into[next_node] = neighbor.edge_id;
                node_queue.emplace(next_node, dist_to_next);
            }
        }
    }
};

/// @brief Computes the shortest path from source to sink using Dijkstra's algorithm.
/// @param graph Current graph.
/// @param source_id ID of source node for path.
/// @param sink_id  ID of final node for path.
/// @param ignored_edges IDs of edges to ignore when pathfinding.
/// @param ignored_nodes IDs of nodes to ignore when pathfinding.
/// @return Path from source to sink.
Path<> Dijkstra(
    const mg_graph::GraphView<> &graph, std::uint64_t source_id, std::uint64_t sink_id,
    const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
    const CheckAbortFunc& check_abort
);

} // namespace shortest_paths