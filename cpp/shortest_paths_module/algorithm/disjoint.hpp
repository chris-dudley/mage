#pragma once

#include <limits>
#include <algorithm>

#include "shortest_path.hpp"
#include "dijkstra.hpp"

namespace shortest_paths {

/// @brief Pathfinder that searches for disjoint shortest paths between a source and target.
/// @tparam TSize Type used for node and edge IDs.
template <typename TSize = std::uint64_t>
class DisjointPathfinder {
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
    /// @brief A vector of edge scores.
    using ScoreVec = std::vector<double>;

private:
    using DijkstraPF = DijkstraPathfinder<TSize>;

    // Combines an Edge Id and score for partial disjoint paths, with ordering being done only on
    // the score.
    struct EdgeAndScore {
        TSize edge_id;
        double score;

        std::partial_ordering operator<=>(const EdgeAndScore& rhs) const {
            return score <=> rhs.score;
        }
    };

public:
    /// @brief Search for the disjoint K-shortest paths.
    /// @param graph The Graph to search.
    /// @param source_id ID of source node for pathfinding.
    /// @param target_id ID of target node for pathfinding.
    /// @param K Number of shortest paths to search for. If 0, find all disjoint paths.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @return A vector containing the shortest paths in order of increasing total cost.
    PathVec search(
        const GraphViewType& graph, TSize source_id, TSize target_id, size_t K,
        const CheckAbortFunc& check_abort = CheckAbortNoop
    ) {
        return search(graph, source_id, target_id, K, {}, {}, check_abort);
    }

    /// @brief Search for the disjoint K-shortest paths.
    /// @param graph The Graph to search.
    /// @param source_id ID of source node for pathfinding.
    /// @param target_id ID of target node for pathfinding.
    /// @param K Number of shortest paths to search for. If 0, find all disjoint paths.
    /// @param ignored_edges IDs of edges to ignore during pathfinding.
    /// @param ignored_nodes IDs of nodes to ignore during pathfinding.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @return A vector containing the shortest paths in order of increasing total cost.
    PathVec search(
        const GraphViewType& graph, TSize source_id, TSize target_id, size_t K,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort = CheckAbortNoop
    ) {
        // If K = 0, find all disjoint paths.
        if (K == 0) {
            K = std::numeric_limits<size_t>::max();
        }

        PathVec result;
        
        DijkstraPF pathfinder;
        EdgeIdSet cur_ignored_edges(ignored_edges);

        for (size_t k = 0; k < K; k++) {
            check_abort();
            pathfinder.search(graph, source_id, target_id, cur_ignored_edges, ignored_nodes, check_abort);
            if (!pathfinder.has_path_to(target_id)) {
                break;
            }
            auto path = pathfinder.path_to(target_id);
            result.push_back(path);

            // Ignore edges in this path for future searches.
            for (auto edge_id : path.edges) {
                cur_ignored_edges.emplace(edge_id);
            }
        }

        return result;
    }

    /// @brief Constructs a set of partial disjoint paths, by using a heuristic based on edge scores
    /// to remove edges from previous shortest paths from further pathfinding attempts.
    ///
    /// Each round, the edges in the current shortest path will be ranked by their score, and the
    /// first `edges_removed_per_round` edges will be excluded from future rounds.
    /// @param graph The Graph to search.
    /// @param source_id ID of source node for pathfinding.
    /// @param target_id ID of target node for pathfinding.
    /// @param K Number of shortest paths to search for. If 0, find all disjoint paths.
    /// @param edge_scores Scores used for ranking edges, must be at least as long as the number of edges in the graph.
    /// @param edges_removed_per_round Desired number of edges to cull per round.
    /// @param cull_ascending Whether edges will be culled in order of ascending or descending score.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @return A vector containing the shortest paths in order of increasing total cost.
    PathVec partial_disjoint_paths(
        const GraphViewType& graph, TSize source_id, TSize target_id, size_t K,
        const ScoreVec& edge_scores, size_t edges_removed_per_round, bool cull_ascending = true,
        const CheckAbortFunc& check_abort = CheckAbortNoop
    ) {
        EdgeIdSet ignored_edges;
        NodeIdSet ignored_nodes;
        return partial_disjoint_paths(
            graph, source_id, target_id, K,
            ignored_edges, ignored_nodes,
            edge_scores, edges_removed_per_round, cull_ascending,
            check_abort
        );
    }

    /// @brief Constructs a set of partial disjoint paths, by using a heuristic based on edge scores
    /// to remove edges from previous shortest paths from further pathfinding attempts.
    ///
    /// Each round, the edges in the current shortest path will be ranked by their score, and the
    /// first `edges_removed_per_round` edges will be excluded from future rounds.
    /// @param graph The Graph to search.
    /// @param source_id ID of source node for pathfinding.
    /// @param target_id ID of target node for pathfinding.
    /// @param K Number of shortest paths to search for. If 0, find all disjoint paths.
    /// @param ignored_edges IDs of edges to ignore during pathfinding.
    /// @param ignored_nodes IDs of nodes to ignore during pathfinding.
    /// @param edge_scores Scores used for ranking edges, must be at least as long as the number of edges in the graph.
    /// @param edges_removed_per_round Desired number of edges to cull per round.
    /// @param cull_ascending Whether edges will be culled in order of ascending or descending score.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @return A vector containing the shortest paths in order of increasing total cost.
    PathVec partial_disjoint_paths(
        const GraphViewType& graph, TSize source_id, TSize target_id, size_t K,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const ScoreVec& edge_scores, size_t edges_removed_per_round, bool cull_ascending = true,
        const CheckAbortFunc& check_abort = CheckAbortNoop
    ) {
        // If K = 0, find all partial disjoint paths.
        if (K == 0) {
            K = std::numeric_limits<size_t>::max();
        }

        if (edge_scores.size() < graph.Edges().size()) {
            throw std::invalid_argument("Scores not provided for all edges");
        }

        PathVec result;
        DijkstraPF pathfinder;
        EdgeIdSet cur_ignored_edges(ignored_edges);

        for (size_t k = 0; k < K; k++) {
            check_abort();
            pathfinder.search(graph, source_id, target_id, cur_ignored_edges, ignored_nodes, check_abort);
            if (!pathfinder.has_path_to(target_id)) {
                break;
            }
            auto path = pathfinder.path_to(target_id);
            result.push_back(path);

            // Determine edges to remove
            std::vector<EdgeAndScore> scored_edges;
            for (auto edge_id : path.edges) {
                scored_edges.emplace_back(edge_id, edge_scores[edge_id]);
            }

            if (cull_ascending) {
                std::sort(scored_edges.begin(), scored_edges.end(), std::less<EdgeAndScore>());
            } else {
                std::sort(scored_edges.begin(), scored_edges.end(), std::greater<EdgeAndScore>());
            }

            // Remove up to `edges_removed_per_round` edges from further pathfinding attempts.
            const size_t num_to_remove = std::min(edges_removed_per_round, scored_edges.size());
            for (size_t i = 0; i < num_to_remove; i++) {
                cur_ignored_edges.emplace(scored_edges[i].edge_id);
            }
        }

        return result;
    }
};

/// @brief Search for the disjoint K-shortest paths.
/// @tparam TSize Type used for node and edge IDs.
/// @param graph The Graph to search.
/// @param source_id ID of source node for pathfinding.
/// @param target_id ID of target node for pathfinding.
/// @param K Number of shortest paths to search for. If 0, find all disjoint paths.
/// @param ignored_edges IDs of edges to ignore during pathfinding.
/// @param ignored_nodes IDs of nodes to ignore during pathfinding.
/// @param check_abort Function that should throw an exception if execution should be aborted.
/// @return A vector containing the shortest paths in order of increasing total cost.
template <typename TSize = std::uint64_t>
std::vector<Path<TSize>> DisjointKShortestPaths(
    const mg_graph::GraphView<TSize>& graph, TSize source_id, TSize target_id, size_t K,
    const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
    const CheckAbortFunc& check_abort = CheckAbortNoop
) {
    DisjointPathfinder<TSize> pathfinder;
    return pathfinder.search(graph, source_id, target_id, K, ignored_edges, ignored_nodes, check_abort);
}

/// @brief Search for the disjoint K-shortest paths.
/// @tparam TSize Type used for node and edge IDs.
/// @param graph The Graph to search.
/// @param source_id ID of source node for pathfinding.
/// @param target_id ID of target node for pathfinding.
/// @param K Number of shortest paths to search for. If 0, find all disjoint paths.
/// @param check_abort Function that should throw an exception if execution should be aborted.
/// @return A vector containing the shortest paths in order of increasing total cost.
template <typename TSize = std::uint64_t>
std::vector<Path<TSize>> DisjointKShortestPaths(
    const mg_graph::GraphView<TSize>& graph, TSize source_id, TSize target_id, size_t K,
    const CheckAbortFunc& check_abort = CheckAbortNoop
) {
    DisjointPathfinder<TSize> pathfinder;
    return pathfinder.search(graph, source_id, target_id, K, check_abort);
}

/// @brief Constructs a set of partial disjoint paths, by using a heuristic based on edge scores
/// to remove edges from previous shortest paths from further pathfinding attempts.
///
/// Each round, the edges in the current shortest path will be ranked by their score, and the
/// first `edges_removed_per_round` edges will be excluded from future rounds.
/// @param graph The Graph to search.
/// @param source_id ID of source node for pathfinding.
/// @param target_id ID of target node for pathfinding.
/// @param K Number of shortest paths to search for. If 0, find all disjoint paths.
/// @param edge_scores Scores used for ranking edges, must be at least as long as the number of edges in the graph.
/// @param edges_removed_per_round Desired number of edges to cull per round.
/// @param cull_ascending Whether edges will be culled in order of ascending or descending score.
/// @param check_abort Function that should throw an exception if execution should be aborted.
/// @return A vector containing the shortest paths in order of increasing total cost.
template <typename TSize = std::uint64_t>
std::vector<Path<TSize>> PartialDisjointKShortestPaths(
    const mg_graph::GraphView<TSize>& graph, TSize source_id, TSize target_id, size_t K,
    const std::vector<double>& edge_scores, size_t edges_removed_per_round, bool cull_ascending = true,
    const CheckAbortFunc& check_abort = CheckAbortNoop
) {
    DisjointPathfinder<TSize> pathfinder;
    return pathfinder.partial_disjoint_paths(
        graph, source_id, target_id, K,
        edge_scores, edges_removed_per_round, cull_ascending,
        check_abort
    );
}

/// @brief Constructs a set of partial disjoint paths, by using a heuristic based on edge scores
/// to remove edges from previous shortest paths from further pathfinding attempts.
///
/// Each round, the edges in the current shortest path will be ranked by their score, and the
/// first `edges_removed_per_round` edges will be excluded from future rounds.
/// @param graph The Graph to search.
/// @param source_id ID of source node for pathfinding.
/// @param target_id ID of target node for pathfinding.
/// @param K Number of shortest paths to search for. If 0, find all disjoint paths.
/// @param ignored_edges IDs of edges to ignore during pathfinding.
/// @param ignored_nodes IDs of nodes to ignore during pathfinding.
/// @param edge_scores Scores used for ranking edges, must be at least as long as the number of edges in the graph.
/// @param edges_removed_per_round Desired number of edges to cull per round.
/// @param cull_ascending Whether edges will be culled in order of ascending or descending score.
/// @param check_abort Function that should throw an exception if execution should be aborted.
/// @return A vector containing the shortest paths in order of increasing total cost.
template <typename TSize = std::uint64_t>
std::vector<Path<TSize>> PartialDisjointKShortestPaths(
    const mg_graph::GraphView<TSize>& graph, TSize source_id, TSize target_id, size_t K,
    const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
    const std::vector<double>& edge_scores, size_t edges_removed_per_round, bool cull_ascending = true,
    const CheckAbortFunc& check_abort = CheckAbortNoop
) {
    DisjointPathfinder<TSize> pathfinder;
    return pathfinder.partial_disjoint_paths(
        graph, source_id, target_id, K,
        ignored_edges, ignored_nodes,
        edge_scores, edges_removed_per_round, cull_ascending,
        check_abort
    );
}

} // namespace shortest_paths