#pragma once

#include <limits>
#include <vector>
#include <queue>
#include <cstdint>
#include <optional>
#include <ranges>
#include <unordered_set>

#include <mg_graph.hpp>

#include "path.hpp"
#include "cycles.hpp"

namespace shortest_paths {

/// @brief Class implementing the Bellman-Ford Single Source Shortest Paths pathfinding algorithm.
/// @tparam TSize Type used for node and edge IDs.
template<typename TSize = std::uint64_t>
class BellmanFordPathfinder {
public:
    /// @brief Type of graph view this object expects.
    using GraphViewType = mg_graph::GraphView<TSize>;
    using Self = BellmanFordPathfinder<TSize>;
    using EdgeIdSet = std::unordered_set<TSize>;
    using NodeIdSet = std::unordered_set<TSize>;
    using NodeIdVec = std::vector<TSize>;
    using EdgeIdVec = std::vector<TSize>;

    using GraphType = mg_graph::Graph<TSize>;

    struct Edge {
        TSize id;
        TSize from;
        TSize to;
        double weight;
    };

    using PredecessorVec = std::vector<std::optional<Edge>>;

    static constexpr const double EPSILON = 1e-14;
    static constexpr const double POS_INF = std::numeric_limits<double>::infinity();

private:
    /// @brief Number of verticies in graph being considered.
    size_t num_vertex;
    /// @brief ID of source node for paths.
    TSize source_id;
    /// @brief `dist_to[v]` = distance of shortest path from source to `v`
    std::vector<double> dist_to;
    /// @brief `edge_into[v]` = edge leading into `v` on shortest path from source
    std::vector<std::optional<Edge>> edge_into;
    /// @brief `on_queue[v]` = is vertex `v` currently on the queue?
    std::vector<bool> on_queue;
    /// @brief Queue of verticies to relax.
    std::deque<TSize> queue;
    /// @brief Number of edge relaxations that have taken place.
    size_t relaxations;
    /// @brief Holds a negative cycle if one is found.
    std::optional<Path<TSize>> cycle;

public:
    BellmanFordPathfinder() = default;
    BellmanFordPathfinder(const Self&) = default;
    BellmanFordPathfinder(Self&&) = default;
    Self& operator=(const Self&) = default;
    Self& operator=(Self&&) = default;

    /// @brief Searches for the shortests paths on `graph` from `source` to all other reachable nodes.
    ///     If a negative cycle is detected, that will be stored.
    /// @param graph The graph to search.
    /// @param source The ID of the source node.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @throws std::invalid_argument if the source node is not in the graph.
    BellmanFordPathfinder(const GraphViewType& graph, TSize source, const CheckAbortFunc& check_abort = CheckAbortNoop):
        num_vertex(graph.Nodes().size()),
        source_id(source),
        dist_to(graph.Nodes().size(), POS_INF),
        edge_into(graph.Nodes().size(), std::nullopt),
        on_queue(graph.Nodes().size(), false),
        queue(),
        relaxations(0),
        cycle(std::nullopt)
    {
        if (num_vertex == 0) {
            // Empty graph, just return
            return;
        }
        if (source >= num_vertex) {
            throw std::invalid_argument("source node not in graph");
        }
        EdgeIdSet empty_edges;
        NodeIdSet empty_nodes;
        do_search(graph, source, empty_edges, empty_nodes, check_abort);
    }

    /// @brief Searches for the shortests paths on `graph` from `source` to all other reachable nodes.
    /// 
    /// If a negative cycle is detected, that will be stored.
    /// Any edges or nodes in the specified sets will be ignored during pathfinding.
    /// @param graph The graph to search.
    /// @param source The ID of the source node.
    /// @param ignored_edges Set of edge IDs to ignore during pathfinding.
    /// @param ignored_nodes Set of node IDs to ignore during pathfinding.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @throws std::invalid_argument if the source node is not in the graph.
    BellmanFordPathfinder(
        const GraphViewType& graph, TSize source,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort = CheckAbortNoop
    ):
        num_vertex(graph.Nodes().size()),
        source_id(source),
        dist_to(graph.Nodes().size(), POS_INF),
        edge_into(graph.Nodes().size(), std::nullopt),
        on_queue(graph.Nodes().size(), false),
        queue(),
        relaxations(0),
        cycle(std::nullopt)
    {
        if (num_vertex == 0) {
            // Empty graph, just return
            return;
        }
        if (source >= num_vertex) {
            throw std::invalid_argument("source node not in graph");
        }
        do_search(graph, source, ignored_edges, ignored_nodes, check_abort);
    }

    /// @brief Resets the state of the pathfinder and then searches for the shorest paths on `graph` from `source`
    ///     to all other reachable nodes. If a negative cycle is detected, it will be stored.
    /// @param graph The graph to search.
    /// @param source The ID of the source node.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @throws std::invalid_argument if the source node is not in the graph.
    void search(const GraphViewType& graph, TSize source, const CheckAbortFunc& check_abort = CheckAbortNoop) {
        if (source >= graph.Nodes().size()) {
            throw std::invalid_argument("source node not in graph");
        }
        EdgeIdSet empty_edges;
        NodeIdSet empty_nodes;
        reset(graph, source);
        do_search(graph, source, empty_edges, empty_nodes, check_abort);
    }

    /// @brief Resets the state of the pathfinder and then searches for the shorest paths on `graph` from `source`
    ///     to all other reachable nodes.
    ///
    /// If a negative cycle is detected, it will be stored.
    /// Any edges or nodes in the specified sets will be ignored during pathfinding.
    /// @param graph The graph to search.
    /// @param source The ID of the source node.
    /// @param ignored_edges Set of edge IDs to ignore during pathfinding.
    /// @param ignored_nodes Set of node IDs to ignore during pathfinding.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @throws std::invalid_argument if the source node is not in the graph.
    void search(
        const GraphViewType& graph, TSize source,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort = CheckAbortNoop
    ) {
        if (source >= graph.Nodes().size()) {
            throw std::invalid_argument("source node not in graph");
        }
        reset(graph, source);
        do_search(graph, source, ignored_edges, ignored_nodes, check_abort);
    }
    
    /// @brief Reports if a negative cycle was detected.
    /// @return `true` if a negative cycle was found.
    bool has_negative_cycle() const noexcept {
        return cycle.has_value();
    }

    /// @brief Returns a found negative cycle, if any.
    std::optional<Path<TSize>> negative_cycle() const {
        return cycle;
    }

    /// @brief Returns the number of reachable nodes in the graph, including the source.
    size_t num_reachable_nodes() const {
        size_t result = 0;
        for (size_t node_id = 0; node_id < dist_to.size(); node_id++) {
            if (dist_to[node_id] < POS_INF) {
                result++;
            }
        }
        return result;
    }

    /// @brief Returns the set of nodes that are the source or reachable from the source.
    NodeIdVec reachable_nodes() const {
        std::vector<TSize> result;
        for (size_t node_id = 0; node_id < dist_to.size(); node_id++) {
            if (dist_to[node_id] < POS_INF) {
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
            if (!maybe_edge.has_value()) {
                continue;
            }
            result.push_back(maybe_edge.value().id);
        }
        return result;
    }

    /// @brief Returns a vector containing the predecessor edges for each vertex.
    ///
    /// For each vertex `v`, `predecessors[v]` will either be empty or contain information about
    /// the edge leading into `v` in the minimum-weight predecessor tree constructed by Bellman-Ford.
    ///
    /// If `predecessors[v]` has no value, that vertex is not reachable from the source.
    /// @return The predecessors vector.
    const PredecessorVec& predecessors() const {
        return edge_into;
    }

    /// @brief Returns whether a path from the source to the specified vertex exists.
    /// @param vertex The ID of the destination vertex.
    bool has_path_to(TSize vertex) const {
        if (vertex >= num_vertex) {
            return false;
        }
        return dist_to[vertex] < POS_INF;
    }

    /// @brief Returns the lowest-cost path from the source to the specified vertex, if one exists.
    /// @param vertex The destination vertex.
    /// @return The found path, or an empty path if no path exists.
    /// @throws std::logic_error if a negative cycle exists.
    Path<TSize> path_to(TSize vertex) const {
        if (has_negative_cycle()) {
            throw std::logic_error("negative cycle exists");
        }

        Path<TSize> result{source_id};
        if (!has_path_to(vertex)) {
            return result;
        }

        std::vector<Edge> stack;
        for (auto edge = edge_into[vertex]; edge.has_value(); edge = edge_into[edge->from]) {
            stack.push_back(edge.value());
        }

        for (const auto& edge : std::views::reverse(stack)) {
            result.add_edge(edge.id, edge.from, edge.to, dist_to[edge.to]);
        }

        return result;
    }

private:
    void reset(const GraphViewType& graph, TSize source) {
        num_vertex = graph.Nodes().size();
        source_id = source;
        dist_to.assign(num_vertex, POS_INF);
        edge_into.assign(num_vertex, std::nullopt);
        on_queue.assign(num_vertex, false);
        queue.clear();
        relaxations = 0;
        cycle = std::nullopt;
    }

    void do_search(
        const GraphViewType& graph, TSize source,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort
    ) {
        dist_to[source] = 0.0;
        queue.push_back(source);
        on_queue[source] = true;
        while (!queue.empty() && !has_negative_cycle()) {
            check_abort();

            auto next_node = queue.front();
            queue.pop_front();
            on_queue[next_node] = false;
            relax(graph, next_node, ignored_edges, ignored_nodes);
        }
    }

    void relax(const GraphViewType& graph, TSize vertex, const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes) {
        for (const auto& neighbor : graph.OutNeighbours(vertex)) {
            if (ignored_edges.contains(neighbor.edge_id) || ignored_nodes.contains(neighbor.node_id)) {
                continue;
            }
            auto next_node = neighbor.node_id;
            auto edge_weight = graph.IsWeighted() ? graph.GetWeight(neighbor.edge_id) : 1.0;
            auto dist_to_next = dist_to[vertex] + edge_weight;
            if (dist_to[next_node] > (dist_to_next + EPSILON)) {
                dist_to[next_node] = dist_to_next;
                edge_into[next_node] = Edge{neighbor.edge_id, vertex, next_node, edge_weight};
                if (!on_queue[next_node]) {
                    queue.push_back(next_node);
                    on_queue[next_node] = true;
                }
            }

            if ((++relaxations % num_vertex) == 0) {
                if (find_negative_cycle()) {
                    // Found a negative cycle, no more we can do
                    return;
                }
            }
        }
    }

    // Build a predecesor graph and search for a cycle.
    bool find_negative_cycle() {
        GraphType graph;
        for (size_t i = 0; i < num_vertex; i++) graph.CreateNode(i);
        for (size_t vertex = 0; vertex < num_vertex; vertex++) {
            if (edge_into[vertex]) {
                const Edge& edge = edge_into[vertex].value();
                graph.CreateEdge(
                    edge.from, edge.to,
                    mg_graph::GraphType::kDirectedGraph, std::make_optional(edge.id),
                    true, edge.weight
                );
            }
        }

        DigraphCycleFinder finder(graph);
        auto possible_cycle = finder.cycle();
        if (!possible_cycle) {
            return false;
        }
        // Translate the edge IDs back
        auto found_cycle = possible_cycle.value();
        for (size_t i = 0; i < found_cycle.edges.size(); i++) {
            found_cycle.edges[i] = graph.GetMemgraphEdgeId(found_cycle.edges[i]);
        }
        cycle = found_cycle;
        return has_negative_cycle();
    }

};

/// @brief Calculates the shortest path between two nodes using the Bellman-Ford algorithm.
///
/// An empty path will be returned if there is no path from the source to the target, or if a
/// negative cycle exists.
/// @tparam TSize Type of node and edge IDs in the graph.
/// @param graph The current graph.
/// @param source_id ID of source node for pathfinding.
/// @param target_id ID of target node for pathfinding.
/// @param ignored_edges IDs of edges to ignore during pathfinding.
/// @param ignored_nodes IDs of nodes to ignore during pathfinding.
/// @param check_abort Function used to periodically check whether execution should be aborted.
/// @return The path from source to target, which may be empty if none exists.
template<typename TSize = std::uint64_t>
Path<TSize> BellmanFord(
    const mg_graph::GraphView<TSize>& graph, TSize source_id, TSize target_id,
    const std::unordered_set<TSize>& ignored_edges,
    const std::unordered_set<TSize>& ignored_nodes,
    CheckAbortFunc check_abort
) {
    BellmanFordPathfinder<TSize> pathfinder(graph, source_id, ignored_edges, ignored_nodes, check_abort);
    if (pathfinder.has_negative_cycle()) {
        return Path(source_id);
    }
    return pathfinder.path_to(target_id);
}

} // namespace shortest_paths