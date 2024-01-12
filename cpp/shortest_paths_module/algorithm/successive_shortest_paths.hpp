#pragma once

#include <algorithm>
#include <memory>

#include <fmt/core.h>

#include "path.hpp"
#include "dijkstra.hpp"

namespace shortest_paths {

/// Options for conversions of flow amounts over a path.
enum FlowConversion {
    /// No conversions, flow has same units across the path.
    None,
    /// For an edge from A -> B, factor is interpreted as (B / A).
    TargetOverSource,
    /// For an edge from A -> B, factor is interpreted as (A / B).
    SourceOverTarget,
};

template<typename TSize = std::uint64_t>
class SuccessiveShortestPathsPathfinder {
public:
    /// @brief Type of graph view this object expects.
    using GraphViewType = mg_graph::GraphView<TSize>;
    using Self = SuccessiveShortestPathsPathfinder<TSize>;
    /// @brief Represents a set of edges by ID.
    using EdgeIdSet = std::unordered_set<TSize>;
    /// @brief Represents a set of nodes by ID.
    using NodeIdSet = std::unordered_set<TSize>;
    /// @brief A vector of weights or scores for nodes or edges.
    using WeightVec = std::vector<double>;

    /// @brief Path type used by this pathfinder.
    using PathType = Path<TSize>;
    /// @brief A combination of a Path along with the flows along the path nodes.
    using PathAndFlow = std::pair<PathType, WeightVec>;
    /// @brief A vector containing results of the search with paths and their input flow.
    using ResultVec = std::vector<PathAndFlow>;

private:
    using DijkstraPF = DijkstraPathfinder<TSize>;

public:
    /// @brief Search for the lowest-cost, max-flow paths for a given input.
    /// @param graph The Graph to search.
    /// @param source_id ID of source node for pathfinding.
    /// @param target_id ID of target node for pathfinding.
    /// @param flow_in Desired input flow from the source node.
    /// @param capacities The capacities of each edge in the graph.
    /// @param factors The factors of each edge in the graph, used for flow conversion.
    ///     May be empty if no conversion requested.
    /// @param espilson A threshold at which two floating point numbers are considered identical.
    /// @param flow_conversion How flows should be unit converted along the path, if at all.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @return A vector containing the shortest paths in order of increasing total cost and their
    ///     respective input flows.
    ResultVec search(
        const GraphViewType& graph, TSize source_id, TSize target_id, double flow_in,
        const WeightVec& capacities, const WeightVec& factors,
        const double epsilon, const FlowConversion flow_conversion,
        const CheckAbortFunc& check_abort = CheckAbortNoop
    ) {
        EdgeIdSet ignored_edges;
        NodeIdSet ignored_nodes;
        return search(
            graph, source_id, target_id, flow_in,
            capacities, factors,
            epsilon, flow_conversion,
            ignored_edges, ignored_nodes,
            check_abort
        );
    }

    /// @brief Search for the lowest-cost, max-flow paths for a given input.
    /// @param graph The Graph to search.
    /// @param source_id ID of source node for pathfinding.
    /// @param target_id ID of target node for pathfinding.
    /// @param flow_in Desired input flow from the source node.
    /// @param capacities The capacities of each edge in the graph.
    /// @param factors The factors of each edge in the graph, used for flow conversion.
    ///     May be empty if no conversion requested.
    /// @param espilson A threshold at which two floating point numbers are considered identical.
    /// @param flow_conversion How flows should be unit converted along the path, if at all.
    /// @param ignored_edges IDs of edges to ignore during pathfinding.
    /// @param ignored_nodes IDs of nodes to ignore during pathfinding.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @return A vector containing the shortest paths in order of increasing total cost and their
    ///     respective input flows.
    ResultVec search(
        const GraphViewType& graph, TSize source_id, TSize target_id, double flow_in,
        const WeightVec& capacities, const WeightVec& factors,
        const double epsilon, const FlowConversion flow_conversion,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort = CheckAbortNoop
    ) {
        ResultVec result;
        const size_t num_edges = graph.Edges().size();

        if (flow_in < epsilon) {
            return result;
        }

        if (capacities.size() < num_edges) {
            throw std::invalid_argument(fmt::format(
                "insufficient capacities: {} given < {} edges in graph",
                capacities.size(), num_edges
            ));
        }

        if (factors.size() < num_edges && flow_conversion != FlowConversion::None) {
            throw std::invalid_argument(fmt::format(
                "insufficient factors: {} given < {} edges in graph",
                factors.size(), num_edges
            ));
        }

        WeightVec weights(num_edges, 1.0);
        if (graph.IsWeighted()) {
            for (TSize i = 0; i < num_edges; i++) {
                weights[i] = graph.GetWeight(i);
            }
        }
        // Initialize current edge capacities from initial capacities
        WeightVec cur_capacities(capacities);

        DijkstraPF pathfinder;
        EdgeIdSet cur_ignored_edges(ignored_edges);

        // Ignore any edges with capacity < epsilon to start with
        for (TSize edge_id = 0; edge_id < cur_capacities.size(); edge_id++) {
            if (cur_capacities[edge_id] < epsilon) {
                cur_ignored_edges.emplace(edge_id);
            }
        }

        WeightVec path_factors;
        WeightVec path_capacities;
        while (flow_in >= epsilon) {
            pathfinder.search_with_weights(graph, source_id, target_id, weights, cur_ignored_edges, ignored_nodes, check_abort);
            auto shortest_path = pathfinder.path_to(target_id);
            if (shortest_path.empty()) {
                // No more shortest paths
                break;
            }

            if (flow_conversion != FlowConversion::None) {
                path_factors.resize(shortest_path.size());
            }
            path_capacities.resize(shortest_path.size());
            for (size_t i = 0; i < shortest_path.size(); i++) {
                if (flow_conversion != FlowConversion::None)
                    path_factors[i] = factors[shortest_path.edges[i]];
                path_capacities[i] = cur_capacities[shortest_path.edges[i]];
            }

            WeightVec flows = calculate_flows(shortest_path, path_factors, path_capacities, flow_in, flow_conversion);
            
            // Input flow for path is first value in result vector
            double path_flow_in = flows[0];

            // Remove flow from current capacities of each edge in the path
            for (size_t i = 0; i < shortest_path.edges.size(); i++) {
                auto edge_id = shortest_path.edges[i];
                cur_capacities[edge_id] -= flows[i];
                // Ignore edges that have reached capacity
                if (cur_capacities[edge_id] < epsilon) {
                    cur_ignored_edges.emplace(edge_id);
                }
            }

            result.emplace_back(std::move(shortest_path), std::move(flows));
            flow_in -= path_flow_in;
        }
        
        return result;
    }

private:
    // Given a set of paths with weights and capacities, and a given desired input flow, determines
    // the actual flow over the path.
    //
    // Will restrict the flow over the path if any edge reaches capacity.
    // Result will be flow at each node.
    // result[i] = input flow to edge[i] and output flow of edge[i-1]
    WeightVec calculate_flows(
        const PathType& path, const WeightVec& path_factors, const WeightVec& path_capacities,
        double flow_in, const FlowConversion conversion
    ) {
        const size_t NUM_EDGES = path.edges.size();
        const size_t NUM_NODES = NUM_EDGES + 1;

        // Resulting flow at each node of the path.
        // Input flow will be at most `flow_in`.
        // If any intermediate flow would result in any edge capacities being exceeded, the minimum exceeded
        // capacity will be the new flow.
        WeightVec result(NUM_NODES, 0.0);

        // Set of conversion factors converting a unit at each node back to the path source unit.
        // e.g. for a path A -> B -> C -> D, will contain [1, B->A, C->A, D->A]
        WeightVec conv_to_source(NUM_NODES, 1.0);
        switch (conversion) {
        case FlowConversion::None:
            // Factors are good as-is (all 1)
            break;
        case FlowConversion::SourceOverTarget:
            // Edge = A -> B, Factor = (A / B),  Flow_B = Flow_A * (1 / Factor_AB)
            // Thus, Given B, Flow_A = Flow_B * Factor_AB
            // For A -> B -> C -> D, Given Flow_C, Flow_A = Flow_C * Factor_AB * Factor_BC
            // Conv_C->A = Conv_B->A * Factor_BC,  Conv_B->A = Conv_A->A * Factor_AB, Conv_A->A = 1 (by definition)
            for (size_t i = 1; i < NUM_NODES; i++) {
                conv_to_source[i] = conv_to_source[i-1] * path_factors[i-1];
            }
            break;
        case FlowConversion::TargetOverSource:
            // Edge = A -> B, Factor = (B / A), Flow_B = Flow_A * Factor_AB
            // Thus, Given B, Flow_A = Flow_B * (1 / Factor_AB)
            // For A -> B -> C -> D, Given Flow_C, Flow_A = Flow_C * (1/Factor_AB) * (1/Factor_BC)
            // Conv_C->A = Conv_B->A * (1/Factor_BC),  Conv_B->A = Conv_A->A * (1/Factor_AB), Conv_A->A = 1 (by definition)
            for (size_t i = 1; i < NUM_NODES; i++) {
                conv_to_source[i] = conv_to_source[i-1] * (1.0 / path_factors[i-1]);
            }
            break;
        }

        // Determine input flow as min(flow_in, path_capacities_in_source_units...)
        double flow = flow_in;
        for (size_t i = 0; i < NUM_EDGES; i++) {
            flow = std::min(flow, path_capacities[i] * conv_to_source[i]);
        }

        for (size_t i = 0; i < NUM_NODES; i++) {
            // Convert path flow to each edge's input unit.
            result[i] = flow / conv_to_source[i];
        }

        return result;
    }
};

} // namespace shortest_paths