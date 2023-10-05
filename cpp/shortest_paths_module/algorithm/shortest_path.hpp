#pragma once

#include <functional>
#include <unordered_set>
#include <vector>

#include <mg_graph.hpp>

#include "path.hpp"

namespace shortest_paths {

// Types used by helper functions

/// @brief A type for representing sets of edge Ids for exclusion during pathfinding.
using EdgeIdSet = std::unordered_set<uint64_t>;
/// @brief A type for representing sets of node Ids for exclusion during pathfinding.
using NodeIdSet = std::unordered_set<uint64_t>;

/// @brief Signature of a function used to check if the execution should be aborted.
///     The function is expected to throw an exception if the execution should be aborted,
///     and do nothing otherwise.
using CheckAbortFunc = std::function<void()>;

/// @brief Signature for functions that compute the shortest path between two nodes in a graph.
///    The two set arguments are sets of edges and nodes to ignore, respectively.
using ShortestPathFunc = std::function<
    Path<>(
        const mg_graph::GraphView<> &, std::uint64_t, std::uint64_t,
        const NodeIdSet&, const EdgeIdSet&,
        const CheckAbortFunc&
    )
>;

/// @brief No-op function for passing in to pathfinders when no checking is required.
void CheckAbortNoop();

} // namespace shortest_paths