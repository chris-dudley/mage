#pragma once

#include <cstdint>
#include <optional>
#include <vector>
#include <stdexcept>
#include <ranges>

#include <mg_graph.hpp>

#include "path.hpp"

namespace shortest_paths {

/// @brief Class that looks for cycles in a directed graph.
///
/// If the graph contains at least one cycle, one of the cycles will be recorded and can be
/// retrieved. If the graph is weighted, the weight of the cycle will be returned in the path,
/// otherwise the edge weights will default to 1.
/// @tparam TSize Type of Node and Edge IDs in the graph.
template<typename TSize = std::uint64_t>
class DigraphCycleFinder {
public:
    /// @brief Type of graph view this object expects.
    using GraphViewType = mg_graph::GraphView<TSize>;
    using Self = DigraphCycleFinder<TSize>;

private:
    struct Edge {
        TSize id;
        TSize from;
        TSize to;
        double weight;
    };

    // Has the node been visited?
    std::vector<bool> visited;
    // Previous edge on the path to each vertex.
    std::vector<std::optional<Edge>> edge_into;
    // Whether vertex is currently on the stack for DFS.
    std::vector<bool> on_stack;
    // The found cycle, if any.
    std::optional<Path<TSize>> cycle_;

public:
    /// @brief Default constructor, effectively equivalent to constructing against an
    ///   empty graph.
    DigraphCycleFinder() = default;
    DigraphCycleFinder(const Self&) = default;
    DigraphCycleFinder(Self&&) = default;
    Self& operator=(const Self&) = default;
    Self& operator=(Self&&) = default;

    /// @brief Searches the specified directional graph for cycles.
    /// @param graph The grah to search.
    DigraphCycleFinder(const GraphViewType& graph):
        visited(graph.Nodes().size(), false),
        edge_into(graph.Nodes().size(), std::nullopt),
        on_stack(graph.Nodes().size(), false),
        cycle_(std::nullopt)
    {
        size_t num_vertices = graph.Nodes().size();
        for (TSize cur_vertex = 0; cur_vertex < num_vertices; cur_vertex++) {
            if (!visited[cur_vertex]) {
                dfs(graph, cur_vertex);
            }
        }
    }

    /// @brief Returns true if a negative cycle was found.
    bool has_cycle() const noexcept {
        return cycle_.has_value();
    }

    /// @brief Returns the negative cycle if one was found.
    std::optional<Path<TSize>> cycle() const noexcept {
        return cycle_;
    }

private:
    void dfs(const GraphViewType& graph, TSize vertex) {
        on_stack[vertex] = true;
        visited[vertex] = true;
        
        for (const auto& neighbor : graph.OutNeighbours(vertex)) {
            TSize next_vertex = neighbor.node_id;

            if (cycle_) {
                // Already found a cycle, no reason to continue.
                return;
            }

            if (visited[next_vertex] && on_stack[next_vertex]) {
                // Found a directed cycle, trace it back
                std::vector<Edge> stack;
                double edge_weight = graph.IsWeighted() ? graph.GetWeight(neighbor.edge_id) : 1.0;
                Edge cur_edge{neighbor.edge_id, vertex, next_vertex, edge_weight};
                while (cur_edge.from != next_vertex) {
                    stack.push_back(cur_edge);
                    cur_edge = edge_into[cur_edge.from].value();
                }
                stack.push_back(cur_edge);

                Path<TSize> cycle_path{next_vertex};
                double distance = 0.0;
                for(const auto& edge : std::views::reverse(stack)) {
                    distance += edge.weight;
                    cycle_path.add_edge(edge.id, edge.from, edge.to, distance);
                }

                cycle_ = cycle_path;
                return;
            }

            if (!visited[next_vertex]) {
                double edge_weight = graph.IsWeighted() ? graph.GetWeight(neighbor.edge_id) : 1.0;
                edge_into[next_vertex] = Edge{neighbor.edge_id, vertex, next_vertex, edge_weight};
                dfs(graph, next_vertex);
            }
        }

        on_stack[vertex] = false;
    }
};

} // namespace shortest_paths