#pragma once

#include <stdexcept>
#include <vector>

#include <mg_graph.hpp>

namespace shortest_paths {

/// @brief Contains information about a specific path through the graph.
/// @tparam TSize Integer type used as IDs for nodes and edges.
template <typename TSize = std::uint64_t>
struct Path {
    Path() = delete;
    /// @brief Constructs an empty path from a source with no edges, or costs.
    Path(TSize source): nodes({source}), edges(), costs({0.0}), total_cost(0.0) {};
    /// @brief Contructs a Path from its components by copying. No validation is done on the inputs.
    /// @param nodes List of nodes in the path.
    /// @param edges List of edges taken along the path.
    /// @param costs Cummulative costs at each node along the path.
    /// @param total_cost Total cost of the path.
    Path(const std::vector<TSize>& nodes, const std::vector<TSize>& edges, const std::vector<double>& costs, double total_cost):
        nodes(nodes), edges(edges), costs(costs), total_cost(total_cost)
        {};
    /// @brief Contructs a Path from its components by moving them. No validation is done on the inputs.
    /// @param nodes List of nodes in the path.
    /// @param edges List of edges taken along the path.
    /// @param costs Cummulative costs at each node along the path.
    /// @param total_cost Total cost of the path.
    Path(std::vector<TSize>&& nodes, std::vector<TSize>&& edges, std::vector<double>&& costs, double total_cost):
        nodes(std::move(nodes)), edges(std::move(edges)),
        costs(std::move(costs)), total_cost(total_cost)
        {};

    Path(const Path<TSize>&) = default;
    Path(Path<TSize>&&) = default;
    auto operator<=>(const Path<TSize>&) const = default;
    Path<TSize>& operator=(const Path<TSize>&) = default;
    Path<TSize>& operator=(Path<TSize>&&) = default;

    /// @brief IDs of nodes traversed, starting with the source and ending with the sink.
    std::vector<TSize> nodes;
    /// @brief IDs of edges traversed from node i to i+1
    std::vector<TSize> edges;
    /// @brief Accumulated cost at node i
    std::vector<double> costs;
    /// @brief Total cost of the path
    double total_cost;

    /// @brief Returns true if the path contains no edges.
    /// @return True if the path contains no edges, false otherwise.
    constexpr bool empty() const noexcept;

    /// @brief Returns the number of edges in the path.
    /// @return The number of edges in the path.
    constexpr TSize size() const noexcept;

    /// @brief Adds an edge from a graph view to the path.
    /// @param edge The edge data to add.
    /// @param cummulative_weight The cummulative weight of the path after traveling this edge.
    /// @throws std::invalid_argument if this path is not empty and edge.from does not match the last
    ///     node in this path.
    void add_edge(const mg_graph::Edge<TSize>& edge, double cummulative_weight);

    /// @brief Adds an edge from a graph view to the path.
    /// @param id The ID of the edge being added.
    /// @param from ID of the node this edge leads from. Must be the same as the last node in the path.
    /// @param to ID of the node this edge leads to.
    /// @param cummulative_weight The cummulative weight of the path after traveling this edge.
    /// @throws std::invalid_argument if this path is not empty and edge.from does not match the last
    ///     node in this path.
    void add_edge(TSize id, TSize from, TSize to, double cummulative_weight);

    /// @brief Creates and returns a Path representing the first `length` edges of this path.
    ///     If length is 0, the resulting path will be empty, including the nodes.
    /// @param length The length (in edges) of the path prefix to generate.
    /// @return A Path object representing the prefix.
    /// @throws std::invalid_argument if length is greater than the number of edges in this path.
    Path<TSize> prefix(TSize length) const;

    /// @brief Extends this path by appending the edges and nodes from the other path.
    /// @param other The path with which to extend this path.
    /// @throws std::invalid_argument If this path is not empty and the other path does not
    ///     start with the last node of this path.
    void extend(const Path<TSize>& other);
    
    /// @brief Creates a new path consisting of this paath with `other` appended to the end.
    /// @param other The path to join to the end of this path.
    /// @return The resulting path.
    /// @throws std::invalid_argument If this path is not empty and the other path does not
    ///     start with the last node of this path.
    Path<TSize> join(const Path<TSize>& other) const;

    /// @brief Returns true if `prefix_nodes` is a prefix of this Path's nodes.
    /// @param prefix_nodes A list of node IDs.
    bool has_node_prefix(const std::vector<TSize>& prefix_nodes) const noexcept;
    /// @brief Returns true if the nodes in `prefix` are a prefix of this Path's nodes.
    /// @param prefix The path to check as a prefix.
    bool has_node_prefix(const Path<TSize>& prefix) const noexcept;

    /// @brief Checks if this path contains the edges of `prefix` as a prefix.
    /// @param prefix The path to check as a prefix.
    /// @return true if prefix's edges is a prefix of this path's edges.
    bool has_prefix(const Path<TSize>& prefix) const noexcept;

    /// @brief Checks if this path is a prefix of the other path.
    /// @param other The path to check against.
    /// @return true if this path's edges are a prefix of other's edges.
    bool is_prefix(const Path<TSize>& other) const noexcept;
};

template <typename TSize>
constexpr bool Path<TSize>::empty() const noexcept {
    return edges.empty();
}

template <typename TSize>
constexpr TSize Path<TSize>::size() const noexcept {
    return edges.size();
}

template <typename TSize>
void Path<TSize>::add_edge(const mg_graph::Edge<TSize>& edge, double cummulative_weight) {
    if (nodes.back() != edge.from) {
        throw std::invalid_argument("edge.from does not match last node in path");
    }
    nodes.push_back(edge.to);
    edges.push_back(edge.id);
    costs.push_back(cummulative_weight);
    total_cost = cummulative_weight;
}

template <typename TSize>
void Path<TSize>::add_edge(TSize id, TSize from, TSize to, double cummulative_weight) {
    if (nodes.back() != from) {
        throw std::invalid_argument("edge.from does not match last node in path");
    }
    nodes.push_back(to);
    edges.push_back(id);
    costs.push_back(cummulative_weight);
    total_cost = cummulative_weight;
}


template <typename TSize>
Path<TSize> Path<TSize>::prefix(TSize length) const {
    if (length > size()) {
        throw std::invalid_argument("prefix length greater than path length");
    }

    Path<TSize> result(nodes.front());

    if (length == 0) {
        // Return empty result if length = 0
        return result;
    }

    for (TSize i = 0; i < length; i++) {
        result.edges.push_back(edges[i]);
        result.nodes.push_back(nodes[i+1]);
        result.costs.push_back(costs[i+1]);
    }
    result.total_cost = result.costs[length];

    return result;
}

template <typename TSize>
void Path<TSize>::extend(const Path<TSize>& other) {
    if (empty()) {
        nodes = other.nodes;
        edges = other.edges;
        costs = other.costs;
        total_cost = other.total_cost;
        return;
    }

    if (nodes.back() != other.nodes.front()) {
        throw std::invalid_argument("first node in extending path does not match last node in this path");
    }
    for (size_t i = 1; i < other.nodes.size(); i++) {
        nodes.push_back(other.nodes[i]);
        costs.push_back(other.costs[i] + total_cost);
    }
    edges.insert(edges.end(), other.edges.cbegin(), other.edges.cend());
    total_cost += other.total_cost;
}

template <typename TSize>
Path<TSize> Path<TSize>::join(const Path<TSize>& other) const {
    Path<TSize> result{*this};
    result.extend(other);
    return result;
}

template <typename TSize>
bool Path<TSize>::has_node_prefix(const std::vector<TSize>& prefix_nodes) const noexcept {
    if (prefix_nodes.size() > nodes.size()) {
        return false;
    }

    for (size_t i = 0; i < prefix_nodes.size(); i++) {
        if (nodes[i] != prefix_nodes[i]) {
            return false;
        }
    }

    return true;
}

template <typename TSize>
bool Path<TSize>::has_node_prefix(const Path<TSize>& prefix) const noexcept {
    return has_node_prefix(prefix.nodes);
}

template <typename TSize>
bool Path<TSize>::has_prefix(const Path<TSize>& prefix) const noexcept {
    if (prefix.size() > size() || nodes.front() != prefix.nodes.front()) {
        return false;
    }

    for (size_t i = 0; i < prefix.size(); i++) {
        if (edges[i] != prefix.edges[i]) {
            return false;
        }
    }

    return true;
}

template <typename TSize>
bool Path<TSize>::is_prefix(const Path<TSize>& other) const noexcept {
    return other.has_prefix(*this);
}

} // namespace shortest_paths