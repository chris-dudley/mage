#pragma once

#include <string_view>

namespace shortest_paths {

constexpr const std::string_view kProcedureYens = "yens";
constexpr const std::string_view kProcedureBellmanFord = "bellman_ford";
constexpr const std::string_view kProcedureIterativeBellmanFordTargeted = "iterative_bellman_ford_targeted";
constexpr const std::string_view kProcedureIterativeBellmanFordSssp = "iterative_bellman_ford_sssp";
constexpr const std::string_view kProcedureIterativeBellmanFordSubgraph = "iterative_bellman_ford_subgraph";
constexpr const std::string_view kProcedureJohnsonsSubgraph = "johnsons_subgraph";
constexpr const std::string_view kProcedureJohnsonsSourceSubgraphs = "johnsons_source_subgraphs";
constexpr const std::string_view kProcedureJohnsonsPaths = "johnsons_paths";
constexpr const std::string_view kProcedureJohnsonsKShortest = "johnsons_k_shortest";
constexpr const std::string_view kProcedureJohnsonsDisjointKShortest = "johnsons_disjoint_k_shortest";
constexpr const std::string_view kProcedureDisjointKShortest = "disjoint_k_shortest";

constexpr const std::string_view kArgumentSourceNode = "source_node";
constexpr const std::string_view kArgumentSourceNodes = "source_nodes";
constexpr const std::string_view kArgumentTargetNode = "target_node";
constexpr const std::string_view kArgumentTargetNodes = "target_nodes";
constexpr const std::string_view kArgumentK = "k";
constexpr const std::string_view kArgumentRelationshipWeightProperty = "relationship_weight_property";
constexpr const std::string_view kArgumentDefaultWeight = "default_weight";
constexpr const std::string_view kArgumentRevisitNodes = "revisit_nodes";
constexpr const std::string_view kArgumentRelationshipScoreProperty = "relationship_score_property";
constexpr const std::string_view kArgumentDefaultScore = "default_score";
constexpr const std::string_view kArgumentCullAscending = "cull_ascending";
constexpr const std::string_view kArgumentThreads = "threads";
constexpr const std::string_view kArgumentCullPerRound = "cull_per_round";

constexpr const std::string_view kReturnIndex = "index";
constexpr const std::string_view kReturnSourceNode = "source_node";
constexpr const std::string_view kReturnTargetNode = "target_node";
constexpr const std::string_view kReturnTotalCost = "total_cost";
constexpr const std::string_view kReturnCosts = "costs";
constexpr const std::string_view kReturnPath = "path";
constexpr const std::string_view kReturnNegativeCycle = "negative_cycle";
constexpr const std::string_view kReturnNumEdgesRemoved = "num_edges_removed";
constexpr const std::string_view kReturnRemovedEdges = "removed_edges";
constexpr const std::string_view kReturnNodes = "nodes";
constexpr const std::string_view kReturnEdges = "edges";
constexpr const std::string_view kReturnNodeWeights = "node_weights";
constexpr const std::string_view kReturnEdgeWeights = "edge_weights";

} // namespace shortest_paths