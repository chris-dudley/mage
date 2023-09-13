#pragma once

#include <string_view>
#include <functional>

#include "shortest_path.hpp"

namespace shortest_paths {

/* k_shortest_paths constants */
constexpr const std::string_view kProcedureKShortestPaths = "k_shortest_paths";
constexpr const std::string_view kProcedureSubgraphKShortestPaths = "subgraph_k_shortest_paths";

constexpr const std::string_view kArgumentSubgraphNodes = "subgraph_nodes";
constexpr const std::string_view kArgumentSubgraphEdges = "subgraph_edges";
constexpr const std::string_view kArgumentSourceNode = "source_node";
constexpr const std::string_view kArgumentTargetNode = "target_node";
constexpr const std::string_view kArgumentK = "k";
constexpr const std::string_view kArgumentRelationshipWeightProperty = "relationship_weight_property";
constexpr const std::string_view kArgumentDefaultWeight = "default_weight";

constexpr const std::string_view kReturnIndex = "index";
constexpr const std::string_view kReturnSourceNode = "source_node";
constexpr const std::string_view kReturnTargetNode = "target_node";
constexpr const std::string_view kReturnTotalCost = "total_cost";
constexpr const std::string_view kReturnCosts = "costs";
constexpr const std::string_view kReturnPath = "path";

/// @brief Signature of a function used to check if the execution should be aborted.
///     The function is expected to throw an exception if the execution should be aborted,
///     and do nothing otherwise.
using CheckAbortFunc = std::function<void()>;

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
    CheckAbortFunc check_abort
);

} // namespace shortest_paths