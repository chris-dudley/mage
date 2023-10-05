#pragma once

#include <functional>

#include "shortest_path.hpp"

namespace shortest_paths {

/// @brief Computes the K shortest paths in the graph from source to sink.
/// @param graph Current graph.
/// @param source_id ID of source node for paths.
/// @param sink_id ID of final node for paths.
/// @param K Number of shortest paths to compute.
/// @param shortest_path_func Function used to find the shortest path between two nodes in the graph.
/// @param check_abort Function used to check if execution should be aborted.
/// @return A vector of paths, each contaning edges taken in the path from source to sink.
std::vector<Path<>> KShortestPaths(
    const mg_graph::GraphView<> &graph, std::uint64_t source_id, std::uint64_t sink_id,
    std::uint64_t K, ShortestPathFunc shortest_path_func,
    const CheckAbortFunc& check_abort
);

} // namespace shortest_paths