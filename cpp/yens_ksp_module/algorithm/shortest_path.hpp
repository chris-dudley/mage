#pragma once

#include <functional>
#include <unordered_set>
#include <vector>

#include <mg_graph.hpp>

#include "path.hpp"

namespace yens_alg {

// Types used by helper functions

/// @brief A type for representing sets of edge Ids for exclusion during pathfinding.
using EdgeIdSet = std::unordered_set<uint64_t>;
/// @brief A type for representing sets of node Ids for exclusion during pathfinding.
using NodeIdSet = std::unordered_set<uint64_t>;

/// @brief Signature for functions that compute the shortest path between two nodes in a graph.
///    The two set arguments are sets of edges and nodes to ignore, respectively.
using ShortestPathFunc = std::function<
    Path<>(
        const mg_graph::GraphView<> &, std::uint64_t, std::uint64_t,
        const NodeIdSet&, const EdgeIdSet&
    )
>;

/// @brief Computes the shortest path from source to sink using Dijkstra's algorithm.
/// @param graph Current graph.
/// @param source_id ID of source node for path.
/// @param sink_id  ID of final node for path.
/// @param ignored_edges IDs of edges to ignore when pathfinding.
/// @param ignored_nodes IDs of nodes to ignore when pathfinding.
/// @return Path from source to sink.
Path<> Dijkstra(
    const mg_graph::GraphView<> &graph, std::uint64_t source_id, std::uint64_t sink_id,
    const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes
);

} // namespace yens_alg