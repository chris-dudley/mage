#pragma once

#include <functional>
#include <algorithm>

#include "bellman_ford.hpp"
#include "shortest_path.hpp"

namespace shortest_paths {

template<typename TSize = std::uint64_t>
class IterativeBellmanFordPathfinder {
public:
    /// @brief Type of graph view this object expects.
    using GraphViewType = mg_graph::GraphView<TSize>;
    using Self = BellmanFordPathfinder<TSize>;
    using EdgeIdSet = std::unordered_set<TSize>;
    using EdgeIdVec = std::vector<TSize>;
    using NodeIdSet = std::unordered_set<TSize>;
    using NodeIdVec = std::vector<TSize>;
    using EdgeScoresVec = std::vector<double>;

    using BFPathfinder = BellmanFordPathfinder<TSize>;
    using PredecessorVec = BFPathfinder::PredecessorVec;

private:
    struct ScoredEdge {
        TSize edge_id;
        double score;

        static bool less(const ScoredEdge& lhs, const ScoredEdge& rhs) {
            return lhs.score < rhs.score;
        }

        static bool greater(const ScoredEdge& lhs, const ScoredEdge& rhs) {
            return lhs.score > rhs.score;
        }
    };

    using ScoreComparator = std::function<bool(const ScoredEdge&, const ScoredEdge&)>;

    /// @brief Bellman-Ford pathfinder.
    BFPathfinder pathfinder;

    /// @brief Scores for each edge, used to determine which to remove from the graph to
    /// break cycles.
    EdgeScoresVec stored_edge_scores;

    /// @brief Whether to cull edges in ascending or descending order.
    bool stored_cull_ascending;

    /// @brief IDs of edges removed from the graph to break cycles
    EdgeIdVec _removed_edges;

public:
    IterativeBellmanFordPathfinder(const Self& other):
        pathfinder(other.pathfinder),
        stored_edge_scores(other.stored_edge_scores),
        stored_cull_ascending(other.stored_cull_ascending),
        _removed_edges(other._removed_edges)
    {}
    IterativeBellmanFordPathfinder(Self&& other):
        pathfinder(std::move(other.pathfinder)),
        stored_edge_scores(std::move(other.stored_edge_scores)),
        stored_cull_ascending(other.stored_cull_ascending),
        _removed_edges(std::move(other._removed_edges))
    {}
    Self& operator=(const Self& other) {
        pathfinder = other.pathfinder;
        stored_edge_scores = other.stored_edge_scores;
        stored_cull_ascending = other.stored_cull_ascending;
        _removed_edges = other._removed_edges;
    }
    Self& operator=(Self&& other) {
        pathfinder = std::move(other.pathfinder);
        stored_edge_scores = std::move(other.stored_edge_scores);
        stored_cull_ascending = other.stored_cull_ascending;
        _removed_edges = std::move(other._removed_edges);
    }

    /// @brief Initialize the pathfinder with no scores and the culling order set to ascending.
    IterativeBellmanFordPathfinder(): 
        pathfinder(), stored_edge_scores(), stored_cull_ascending(true),
        _removed_edges()
    {}

    /// @brief Initialize the pathfinder with a set of stored edge scores and culling order.
    /// @param edge_scores The edge scores to store.
    /// @param cull_ascending The culling order to use when removing edges in negative cycles.
    IterativeBellmanFordPathfinder(const EdgeScoresVec& edge_scores, bool cull_ascending = true):
        pathfinder(),
        stored_edge_scores(edge_scores),
        stored_cull_ascending(cull_ascending),
        _removed_edges()
    { }

    /// @brief Sets the stored edge scores to the specified scores.
    /// @param scores New edge scores.
    void set_scores(const EdgeScoresVec& scores) {
        stored_edge_scores = scores;
    }

    /// @brief Sets the stored edge score comparator to cull edges in ascending or descending
    /// order based on score.
    /// @param ascending If true, cull edges in ascending order of score. If false, in descending order.
    void cull_ascending(bool ascending) {
        stored_cull_ascending = ascending;
    }

    /// @brief Returns the number of edges that were excluded during pathfinding. Roughly equivalent to
    /// the number of cycles that were broken.
    size_t num_edges_removed() const {
        return _removed_edges.size();
    }

    /// @brief Returns the IDs of edges removed from the graph during pathfinding to break cycles.
    const EdgeIdVec& removed_edges() const {
        return _removed_edges;
    }

    /// @brief Returns the number of reachable nodes in the graph, including the source.
    size_t num_reachable_nodes() const {
        return pathfinder.num_reachable_nodes();
    }

    /// @brief Returns the set of nodes that are the source or reachable from the source.
    NodeIdVec reachable_nodes() const {
        return pathfinder.reachable_nodes();
    }

    /// @brief Returns the set of edges used to traverse the minimum-weighted paths to reachable nodes.
    /// If the graph has a negative cycle, the returned edges are not guaranteed to contain all edges needed
    /// to reach all nodes.
    EdgeIdVec edges_used() const {
        return pathfinder.edges_used();
    }

    /// @brief Returns whether a path from the source to the specified vertex exists.
    /// @param vertex The ID of the destination vertex.
    bool has_path_to(TSize vertex) const {
        return pathfinder.has_path_to(vertex);
    }

    /// @brief Returns the lowest-cost path from the source to the specified vertex, if one exists.
    /// @param vertex The destination vertex.
    /// @return The found path, or an empty path if no path exists.
    /// @throws std::logic_error if a negative cycle exists.
    Path<TSize> path_to(TSize vertex) const {
        return pathfinder.path_to(vertex);
    }

    /// @brief Returns a vector containing the predecessor edges for each vertex.
    ///
    /// For each vertex `v`, `predecessors[v]` will either be empty or contain information about
    /// the edge leading into `v` in the minimum-weight predecessor tree constructed by Bellman-Ford.
    ///
    /// If `predecessors[v]` has no value, that vertex is not reachable from the source.
    /// @return The predecessors vector.
    const PredecessorVec& predecessors() const {
        return pathfinder.predecessors();
    }

    /// @brief Searches the graph for a path from source to target without negative cycles. The previously
    /// stored edge scores and culling order will be used to remove edges if negative cycles are encountered.
    /// @param graph The graph to operate on.
    /// @param source The ID of the source node.
    /// @param target The ID of the target node.
    /// @param initial_ignored_edges Initial set of edges to ignore, if any.
    /// @param ignored_nodes Set of nodes to ignore during pathfinding.
    /// @param check_abort Function used to checked periodically whether execution should be aborted.
    /// @return The path from source to target, which may be empty if none is found.
    /// @throws std::invalid_argument If the previously stored edge scores does not contain scores for all edges in the
    /// graph.
    Path<TSize> search(
        const GraphViewType& graph, TSize source, TSize target,
        const EdgeIdSet& initial_ignored_edges, const NodeIdSet& ignored_nodes,
        CheckAbortFunc check_abort
    ) {
        return search(
            graph, source, target,
            initial_ignored_edges, ignored_nodes,
            stored_edge_scores, stored_cull_ascending,
            check_abort
        );
    }

    /// @brief Searches the graph for a path from source to target without negative cycles.
    ///
    /// Will use the given edge scores and culling order to remove edges if negative cycles are encountered.
    /// @param graph The graph to operate on.
    /// @param source The ID of the source node.
    /// @param target The ID of the target node.
    /// @param initial_ignored_edges Initial set of edges to ignore, if any.
    /// @param ignored_nodes Set of nodes to ignore during pathfinding.
    /// @param edge_scores The scores for each edge to be used when removing edges in negative cycles.
    /// @param cull_ascending Whether culling is performed in ascending or descending order of score.
    /// @param check_abort Function used to checked periodically whether execution should be aborted.
    /// @return The path from source to target, which may be empty if none is found.
    /// @throws std::invalid_argument If the edge scores does not contain scores for all edges in the
    /// graph.
    Path<TSize> search(
        const GraphViewType& graph, TSize source, TSize target,
        const EdgeIdSet& initial_ignored_edges, const NodeIdSet& ignored_nodes,
        const EdgeScoresVec& edge_scores, bool cull_ascending,
        CheckAbortFunc check_abort
    ) {
        Path<TSize> result(source);
        if (source == target) {
            return result;
        }

        auto num_edges = graph.Edges().size();
        if (num_edges > edge_scores.size()) {
            throw std::invalid_argument("Scores not available for all edges in graph");
        }

        do_search(
            graph, source, target,
            initial_ignored_edges, ignored_nodes,
            edge_scores, cull_ascending, check_abort
        );

        if (pathfinder.has_path_to(target)) {
            result = pathfinder.path_to(target);
        }
        return result;
    }

    /// @brief Constructs paths from the source node to all reachable nodes, while removing cycles.
    ///
    /// The previously stored edge scores and culling order will be used to determine which edges to
    /// remove.
    /// @param graph The graph to operate on.
    /// @param source The ID of the source node.
    /// @param initial_ignored_edges Initial set of edges to ignore, if any.
    /// @param ignored_nodes Set of nodes to ignore during pathfinding.
    /// @param check_abort Function used to checked periodically whether execution should be aborted.
    /// @throws std::invalid_argument If the stored edge scores do not contain scores for all edges in the
    /// graph.
    void remove_cycles(
        const GraphViewType& graph, TSize source,
        const EdgeIdSet& initial_ignored_edges, const NodeIdSet& ignored_nodes,
        CheckAbortFunc check_abort
    ) {
        remove_cycles(
            graph, source,
            initial_ignored_edges, ignored_nodes,
            stored_edge_scores, stored_cull_ascending,
            check_abort
        );
    }

    /// @brief Constructs paths from the source node to all reachable nodes, while removing cycles.
    ///
    /// Will use the given edge scores and culling order to remove edges if negative cycles are encountered.
    /// @param graph The graph to operate on.
    /// @param source The ID of the source node.
    /// @param initial_ignored_edges Initial set of edges to ignore, if any.
    /// @param ignored_nodes Set of nodes to ignore during pathfinding.
    /// @param edge_scores The scores for each edge to be used when removing edges in negative cycles.
    /// @param cull_ascending Whether culling is performed in ascending or descending order of score.
    /// @param check_abort Function used to checked periodically whether execution should be aborted.
    /// @throws std::invalid_argument If the stored edge scores do not contain scores for all edges in the
    /// graph.
    void remove_cycles(
        const GraphViewType& graph, TSize source,
        const EdgeIdSet& initial_ignored_edges, const NodeIdSet& ignored_nodes,
        const EdgeScoresVec& edge_scores, bool cull_ascending,
        CheckAbortFunc check_abort
    ) {
        auto num_edges = graph.Edges().size();
        if (num_edges > edge_scores.size()) {
            throw std::invalid_argument("Scores not available for all edges in graph");
        }
        do_search_no_target(
            graph, source,
            initial_ignored_edges, ignored_nodes,
            edge_scores, cull_ascending,
            check_abort
        );
    }

private:
    void do_search(
        const GraphViewType& graph, TSize source, TSize target,
        const EdgeIdSet& initial_ignored_edges, const NodeIdSet& ignored_nodes,
        const EdgeScoresVec& edge_scores, bool cull_ascending,
        CheckAbortFunc check_abort
    ) {
        ScoreComparator comparator = cull_ascending ? ScoredEdge::greater : ScoredEdge::less;
        EdgeIdSet ignored_edges = initial_ignored_edges;
        _removed_edges.clear();

        pathfinder.search(graph, source, ignored_edges, ignored_nodes, check_abort);
        if (!pathfinder.has_path_to(target)) {
            // No path exists to target
            return;
        }

        while (pathfinder.has_negative_cycle()) {
            check_abort();

            auto cycle = pathfinder.negative_cycle().value();

            // Order the edges in the cycle by score and try to remove
            std::vector<ScoredEdge> cycle_edges;
            for (auto cycle_edge_id : cycle.edges) {
                cycle_edges.emplace_back(cycle_edge_id, edge_scores[cycle_edge_id]);
            }
            std::sort(cycle_edges.begin(), cycle_edges.end(), comparator);

            bool made_progress = false;
            while (!cycle_edges.empty()) {
                auto current_edge = cycle_edges.back();
                cycle_edges.pop_back();

                ignored_edges.emplace(current_edge.edge_id);
                pathfinder.search(graph, source, ignored_edges, ignored_nodes, check_abort);

                if (pathfinder.has_negative_cycle() || pathfinder.has_path_to(target)) {
                    // We either found a new negative cycle, or removing that edge found us a path to the target.
                    made_progress = true;
                    _removed_edges.push_back(current_edge.edge_id);
                    break;
                }

                // Removing that edge broke the path from source to target, have to try again
                ignored_edges.erase(current_edge.edge_id);
            }

            if (!made_progress) {
                throw std::logic_error("Unable to break cycle without breaking path to target");
            }
        }
    }

    void do_search_no_target(
        const GraphViewType& graph, TSize source,
        const EdgeIdSet& initial_ignored_edges, const NodeIdSet& ignored_nodes,
        const EdgeScoresVec& edge_scores, bool cull_ascending,
        CheckAbortFunc check_abort
    ) {
        // Returns true if lhs is worse than rhs
        ScoreComparator worse_than = cull_ascending ? ScoredEdge::greater : ScoredEdge::less;

        EdgeIdSet ignored_edges = initial_ignored_edges;
        _removed_edges.clear();

        pathfinder.search(graph, source, ignored_edges, ignored_nodes, check_abort);

        while (pathfinder.has_negative_cycle()) {
            check_abort();

            auto cycle = pathfinder.negative_cycle().value();

            ScoredEdge current_edge(cycle.edges[0], edge_scores[cycle.edges[0]]);
            for (size_t i = 1; i < cycle.edges.size(); i++) {
                ScoredEdge next_edge(cycle.edges[i], edge_scores[cycle.edges[i]]);
                if (worse_than(current_edge, next_edge)) {
                    current_edge = next_edge;
                }
            }

            ignored_edges.emplace(current_edge.edge_id);
            _removed_edges.push_back(current_edge.edge_id);

            pathfinder.search(graph, source, ignored_edges, ignored_nodes, check_abort);
        }
    }
};

/// @brief Calculates the shortest path between two nodes using the modified iterative Bellman-Ford algorithm.
///
/// The scores will default to the edge weights, with culling in descending order.
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
Path<TSize> IterativeBellmanFord(
    const mg_graph::GraphView<TSize>& graph, TSize source_id, TSize target_id,
    const std::unordered_set<TSize>& ignored_edges,
    const std::unordered_set<TSize>& ignored_nodes,
    CheckAbortFunc check_abort
) {
    size_t num_edges = graph.Edges().size();
    std::vector<double> scores(num_edges, 1.0);
    if (graph.IsWeighted()) {
        for (TSize edge_id = 0; edge_id < num_edges; edge_id++) {
            scores[edge_id] = graph.GetWeight(edge_id);
        }
    }
    IterativeBellmanFordPathfinder<TSize> pathfinder(scores, false);
    Path<TSize> path = pathfinder.search(graph, source_id, target_id, ignored_edges, ignored_nodes, check_abort);
    return path;
}

/// @brief Calculates the shortest path between two nodes using the modified iterative Bellman-Ford algorithm.
///
/// An empty path will be returned if there is no path from the source to the target, or if a
/// negative cycle exists.
/// @tparam TSize Type of node and edge IDs in the graph.
/// @param graph The current graph.
/// @param source_id ID of source node for pathfinding.
/// @param target_id ID of target node for pathfinding.
/// @param ignored_edges IDs of edges to ignore during pathfinding.
/// @param ignored_nodes IDs of nodes to ignore during pathfinding.
/// @param edge_scores Scores for each edge, used to determine the order to attempt to cull edges.
/// @param cull_ascending Whether edges should be culled in order of ascending or descending score.
/// @param check_abort Function used to periodically check whether execution should be aborted.
/// @return The path from source to target, which may be empty if none exists.
template<typename TSize = std::uint64_t>
Path<TSize> IterativeBellmanFord(
    const mg_graph::GraphView<TSize>& graph, TSize source_id, TSize target_id,
    const std::unordered_set<TSize>& ignored_edges,
    const std::unordered_set<TSize>& ignored_nodes,
    const std::vector<double>& edge_scores,
    bool cull_ascending,
    CheckAbortFunc check_abort
) {
    size_t num_edges = graph.Edges().size();
    if (num_edges > edge_scores.size()) {
        throw std::invalid_argument("Scores not provided for all edges");
    }
    IterativeBellmanFordPathfinder<TSize> pathfinder(edge_scores, cull_ascending);
    Path<TSize> path = pathfinder.search(graph, source_id, target_id, ignored_edges, ignored_nodes, check_abort);
    return path;
}

} // namespace shortest_paths