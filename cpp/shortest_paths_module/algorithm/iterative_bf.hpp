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
    using NodeIdSet = std::unordered_set<TSize>;
    using EdgeScoresVec = std::vector<double>;

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
    using BFPathfinder = BellmanFordPathfinder<TSize>;

    /// @brief Bellman-Ford pathfinder.
    BFPathfinder pathfinder;

    /// @brief Scores for each edge, used to determine which to remove from the graph to
    /// break cycles.
    EdgeScoresVec stored_edge_scores;

    /// @brief Whether to cull edges in ascending or descending order.
    bool stored_cull_ascending;

public:
    IterativeBellmanFordPathfinder(const Self& other):
        pathfinder(other.pathfinder),
        stored_edge_scores(other.stored_edge_scores),
        stored_cull_ascending(other.stored_cull_ascending)
    {}
    IterativeBellmanFordPathfinder(Self&& other):
        pathfinder(std::move(other.pathfinder)),
        stored_edge_scores(std::move(other.stored_edge_scores)),
        stored_cull_ascending(other.stored_cull_ascending)
    {}
    Self& operator=(const Self& other) {
        pathfinder = other.pathfinder;
        stored_edge_scores = other.stored_edge_scores;
        stored_cull_ascending = other.stored_cull_ascending;
    }
    Self& operator=(Self&& other) {
        pathfinder = std::move(other.pathfinder);
        stored_edge_scores = std::move(other.stored_edge_scores);
        stored_cull_ascending = other.stored_cull_ascending;
    }

    /// @brief Initialize the pathfinder with no scores and the culling order set to ascending.
    IterativeBellmanFordPathfinder(): pathfinder(), stored_edge_scores(), stored_cull_ascending(true) {}

    /// @brief Initialize the pathfinder with a set of stored edge scores and culling order.
    /// @param edge_scores The edge scores to store.
    /// @param cull_ascending The culling order to use when removing edges in negative cycles.
    IterativeBellmanFordPathfinder(const EdgeScoresVec& edge_scores, bool cull_ascending = true):
        pathfinder(),
        stored_edge_scores(edge_scores),
        stored_cull_ascending(cull_ascending)
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
        Path<TSize> result(source);
        if (source == target) {
            return result;
        }

        auto num_edges = graph.Edges().size();
        if (num_edges > stored_edge_scores.size()) {
            throw std::invalid_argument("Scores not available for all edges in graph");
        }

        ScoreComparator comparator = stored_cull_ascending ? ScoredEdge::greater : ScoredEdge::less;
        do_search(
            graph, source, target,
            initial_ignored_edges, ignored_nodes,
            stored_edge_scores, comparator, check_abort
        );

        if (pathfinder.has_path_to(target)) {
            result = pathfinder.path_to(target);
        }
        return result;
    }

    /// @brief Searches the graph for a path from source to target without negative cycles. Will use the
    /// given edge scores and culling order to remove edges if negative cycles are encountered.
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
        ScoreComparator comparator = cull_ascending ? ScoredEdge::greater : ScoredEdge::less;

        do_search(
            graph, source, target,
            initial_ignored_edges, ignored_nodes,
            edge_scores, comparator, check_abort
        );

        if (pathfinder.has_path_to(target)) {
            result = pathfinder.path_to(target);
        }
        return result;
    }

private:
    void do_search(
        const GraphViewType& graph, TSize source, TSize target,
        const EdgeIdSet& initial_ignored_edges, const NodeIdSet& ignored_nodes,
        const EdgeScoresVec& edge_scores, ScoreComparator comparator,
        CheckAbortFunc check_abort
    ) {
        EdgeIdSet ignored_edges = initial_ignored_edges;

        pathfinder.search(graph, source, ignored_edges, ignored_nodes);

        while (pathfinder.has_negative_cycle()) {
            check_abort();

            auto cycle = pathfinder.negative_cycle().value();

            // Order the edges in the cycle by score and try to remove
            std::vector<ScoredEdge> cycle_edges;
            for (auto cycle_edge_id : cycle.edges) {
                cycle_edges.emplace_back(cycle_edge_id, edge_scores[cycle_edge_id]);
            }
            std::make_heap(cycle_edges.begin(), cycle_edges.end(), comparator);

            while (!cycle_edges.empty()) {
                std::pop_heap(cycle_edges.begin(), cycle_edges.end(), comparator);
                auto current_edge = cycle_edges.back();
                cycle_edges.pop_back();

                ignored_edges.emplace(current_edge.edge_id);
                pathfinder.search(graph, source, ignored_edges, ignored_nodes);

                if (pathfinder.has_negative_cycle() || pathfinder.has_path_to(target)) {
                    // We either found a new negative cycle, or removing that edge found us a path to the target.
                    break;
                }

                // Removing that edge broke the path from source to target, have to try again
                ignored_edges.erase(current_edge.edge_id);
            }
        }
    }
};

} // namespace shortest_paths