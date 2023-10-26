#pragma once

#include <limits>

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

private:
    using DijkstraPF = DijkstraPathfinder<TSize>;

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

} // namespace shortest_paths