#pragma once

#include <string_view>

namespace shortest_paths {

constexpr const std::string_view kProcedureYens = "yens";
constexpr const std::string_view kProcedureYensSubgraph = "yens_subgraph";
constexpr const std::string_view kProcedureBellmanFord = "bellman_ford";
constexpr const std::string_view kProcedureIterativeBellmanFord = "iterative_bellman_ford";

constexpr const std::string_view kArgumentSubgraphNodes = "subgraph_nodes";
constexpr const std::string_view kArgumentSubgraphEdges = "subgraph_edges";
constexpr const std::string_view kArgumentSourceNode = "source_node";
constexpr const std::string_view kArgumentTargetNode = "target_node";
constexpr const std::string_view kArgumentTargetNodes = "target_nodes";
constexpr const std::string_view kArgumentK = "k";
constexpr const std::string_view kArgumentRelationshipWeightProperty = "relationship_weight_property";
constexpr const std::string_view kArgumentDefaultWeight = "default_weight";
constexpr const std::string_view kArgumentRevisitNodes = "revisit_nodes";
constexpr const std::string_view kArgumentRelationshipScoreProperty = "relationship_score_property";
constexpr const std::string_view kArgumentDefaultScore = "default_score";
constexpr const std::string_view kArgumentCullAscending = "cull_ascending";

constexpr const std::string_view kReturnIndex = "index";
constexpr const std::string_view kReturnSourceNode = "source_node";
constexpr const std::string_view kReturnTargetNode = "target_node";
constexpr const std::string_view kReturnTotalCost = "total_cost";
constexpr const std::string_view kReturnCosts = "costs";
constexpr const std::string_view kReturnPath = "path";
constexpr const std::string_view kReturnNegativeCycle = "negative_cycle";
constexpr const std::string_view kReturnEdgesRemoved = "edges_removed";

} // namespace shortest_paths