#pragma once

#include <memory>
#include <thread>
#include <ranges>

#include <fmt/core.h>
#include <omp.h>

#include "path.hpp"
#include "shortest_path.hpp"
#include "bellman_ford.hpp"
#include "iterative_bf.hpp"
#include "dijkstra.hpp"

namespace shortest_paths {

template<typename TSize = std::uint64_t>
class JohnsonsPathfinder {
public:
    /// @brief Type of graph view this object expects.
    using GraphViewType = mg_graph::GraphView<TSize>;
    using Self = JohnsonsPathfinder<TSize>;
    /// @brief Represents a set of edges by ID.
    using EdgeIdSet = std::unordered_set<TSize>;
    /// @brief Represents a set of nodes by ID.
    using NodeIdSet = std::unordered_set<TSize>;
    /// @brief A vector of nodes by ID.
    using NodeIdVec = std::vector<TSize>;
    /// @brief A vector of edges by ID.
    using EdgeIdVec = std::vector<TSize>;
    /// @brief A vector of weights or scores for nodes or edges.
    using WeightVec = std::vector<double>;
    /// @brief A Dijkstra pathfinder, used for finding all paths from a source node.
    using DijkstraPF = DijkstraPathfinder<TSize>;

    /// @brief Infinite distance used for unreachable nodes.
    static constexpr const double POSITIVE_INFINITY = std::numeric_limits<double>::infinity();

private:
    using GraphType = mg_graph::Graph<TSize>;
    using BellmanFordPF = BellmanFordPathfinder<TSize>;
    using IteraveBellmanFordPF = IterativeBellmanFordPathfinder<TSize>;
    using DijkstraUniquePtr = std::unique_ptr<DijkstraPF>;

    // Adjusted nodes weights calculated using Bellman-Ford's algorithm.
    std::vector<double> _node_weights;

    // Adjusted edge weights. Weights are calculated, for an edge `e` between nodes `u` and `v`:
    // `edge_weights[e] = weight(e) + node_weights[u] - node_weights[v]`
    std::vector<double> _edge_weights;

    // Vector of Dijkstra pathfinders, one for each source. Each may or may not be initialized based on
    // what paths are requested.
    std::vector<DijkstraUniquePtr> _pathfinders;

    // Negative cycle, if standard Bellman-Ford pathfinder is used and a negative cycle is found.
    std::optional<Path<TSize>> _neg_cycle;

    // Edges removed when using iterative Bellman-Ford.
    EdgeIdSet _removed_edges;

public:
    JohnsonsPathfinder() = default;
    JohnsonsPathfinder(Self&&) = default;
    Self& operator=(Self&&) = default;

    /// @brief Returns true if the pathfinder encountered a negative cycle and was not instructed
    /// to remove them.
    bool has_negative_cycle() const {
        return _neg_cycle.has_value();
    }

    /// @brief Returns the negative cycle that was encountered and not removed, if any.
    const std::optional<Path<TSize>>& negative_cycle() const {
        return _neg_cycle;
    }

    /// @brief Returns the number of edges removed from the graph in order to break cycles.
    size_t num_edges_removed() const {
        return _removed_edges.size();
    }

    /// @brief Returns the set of edges removed in order to break cycles.
    const EdgeIdSet& removed_edges() const {
        return _removed_edges;
    }

    /// @brief Returns the adjusted node weights calculated using Bellman-Ford.
    const WeightVec& node_weights() const {
        return _node_weights;
    }

    /// @brief Returns the adjusted edge weights, where for an an edge `e` from `u` -> `v`,
    /// `edge_weights[e] = weight(e) + node_weights[u] - node_weights[v]`
    const WeightVec& edge_weights() const {
        return _edge_weights;
    }

    /// @brief Returns whether a path exists from `source` to `target`.
    /// @param source ID of the source node.
    /// @param target ID of the target node.
    /// @return True if a path exists, false if no path exists or pathfinding from source was not
    /// performed.
    bool has_path(TSize source, TSize target) const {
        if (!has_pathfinder(source)) {
            return false;
        }
        return _pathfinders[source]->has_path_to(target);
    }

    /// @brief Returns the path from source to target, if one exists.
    /// @param source ID of source node of path.
    /// @param target ID of target node of path.
    /// @return The path from source to target, which may be empty if no such path exists.
    Path<TSize> get_path(TSize source, TSize target) const {
        if (!has_pathfinder(source)) {
            return Path<TSize>(source);
        }
        return _pathfinders[source]->path_to(target);
    }

    /// @brief Returns the number of nodes that are reachable from the source, including the source.
    /// @param source ID of the source node.
    size_t num_reachable_nodes_from(TSize source) const {
        if (!has_pathfinder(source)) {
            return 0;
        }
        return _pathfinders[source]->num_reachable_nodes();
    }

    /// @brief Returns the set of nodes that are reachable from the specified source
    /// (including the source).
    /// Will be empty if no search from the source was requested.
    /// @param source ID of the source node.
    NodeIdVec reachable_nodes_from(TSize source) const {
        if (!has_pathfinder(source)) {
            return {};
        }
        return _pathfinders[source]->reachable_nodes();
    }

    /// @brief Returns the set of edges used to traverse the minimum-weighted paths to reachable nodes
    /// from the specified source.
    ///
    /// If the graph has a negative cycle, the returned edges are not guaranteed to contain all edges needed
    /// to reach all nodes.
    EdgeIdVec edges_used_from(TSize source) const {
        if (!has_pathfinder(source)) {
            return {};
        }
        return _pathfinders[source]->edges_used();
    }

    /// @brief Returns the distances from the specified source to all nodes in the graph.
    /// @param source ID of the source node.
    /// @return Vector of distances, where `distance[t]` is the shortest distance from `source` to `t`.
    /// @throws std::invalid_argument if pathfinding was not requested for source.
    const std::vector<double>& distances_from(TSize source) const {
        return get_pathfinder(source).distances();
    }

    /// @brief Returns whether the Dijktra pathfinder for the given source has been initialized.
    /// @param source ID of the source node.
    /// @return True if the pathfinder for the source node is initialized and available for queries.
    bool has_pathfinder(TSize source) const {
        if (source >= _pathfinders.size()) {
            return false;
        }
        return bool(_pathfinders[source]);
    }

    /// @brief Returns the underlying Dijkstra pathfinder for the given source node.
    /// @param source ID of the source node.
    /// @return A reference to the pathfinder, which may be queried.
    /// @throws std::invalid_argument if pathfinding was not requested for source.
    const DijkstraPF& get_pathfinder(TSize source) const {
        if (!has_pathfinder(source)) {
            throw std::invalid_argument(fmt::format("Pathfinder not initialized for source node {}", source));
        }
        return *(_pathfinders[source]);
    }

    /// @brief Attempts to compute the all pairs shortest paths for the given graph.
    ///
    /// If a negative cycle is detected, it will be stored and pathfinding will be aborted.
    ///
    /// @param graph The graph to search.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @param threads The maximum number of threads to use during pathfinding. If <= 0, will be set to
    /// the number of processors available.
    void search_all(
        const GraphViewType& graph,
        const CheckAbortFunc& check_abort = CheckAbortNoop,
        int threads = 0
    ) {
        EdgeIdSet ignored_edges;
        NodeIdSet ignored_nodes;
        search_all(graph, ignored_edges, ignored_nodes, check_abort, threads);
    }

    /// @brief Attempts to compute the all pairs shortest paths for the given graph.
    ///
    /// If a negative cycle is detected, it will be stored and pathfinding will be aborted.
    ///
    /// @param graph The graph to search.
    /// @param ignored_edges Set of edge IDs to ignore during pathfinding.
    /// @param ignored_nodes Set of node IDs to ignore during pathfinding.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @param threads The maximum number of threads to use during pathfinding. If <= 0, will be set to
    /// the number of processors available.
    void search_all(
        const GraphViewType& graph,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort = CheckAbortNoop,
        int threads = 0
    ) {
        calculate_node_weights_bf(graph, ignored_edges, ignored_nodes, check_abort);
        if (_neg_cycle) {
            return;
        }
        calculate_edge_weights(graph);

        auto reweighted_graph = copy_graph(graph, _edge_weights);
        if (threads == 1) {
            pathfind_all(reweighted_graph, ignored_edges, ignored_nodes, check_abort);
        } else {
            pathfind_all_omp(reweighted_graph, ignored_edges, ignored_nodes, check_abort, threads);
        }
    }

    /// @brief Attempts to compute the all pairs shortest paths for the given graph.
    ///
    /// If a negative cycle is detected, it will be stored and pathfinding will be aborted.
    ///
    /// @param graph The graph to search.
    /// @param sources A container of source node IDs to pathfind from.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @param threads The maximum number of threads to use during pathfinding. If <= 0, will be set to
    /// the number of processors available.
    /// @tparam TSourceCont A container of source IDs.
    template <typename TSourceCont>
    void search_some(
        const GraphViewType& graph, const TSourceCont& sources,
        const CheckAbortFunc& check_abort = CheckAbortNoop,
        int threads = 0
    ) {
        search_some(graph, sources, {}, {}, check_abort, threads);
    }

    /// @brief Attempts to compute the shortest paths for the given graph from the specified
    /// source nodes.
    ///
    /// If a negative cycle is detected, it will be stored and pathfinding will be aborted.
    ///
    /// @param graph The graph to search.
    /// @param sources A container of source node IDs to pathfind from.
    /// @param ignored_edges Set of edge IDs to ignore during pathfinding.
    /// @param ignored_nodes Set of node IDs to ignore during pathfinding.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @param threads The maximum number of threads to use during pathfinding. If <= 0, will be set to
    /// the number of processors available.
    /// @tparam TSourceCont A container of source IDs.
    template <typename TSourceCont>
    void search_some(
        const GraphViewType& graph, const TSourceCont& sources,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort = CheckAbortNoop,
        int threads = 0
    ) {
        static_assert(std::ranges::range<TSourceCont>);
        static_assert(std::convertible_to< std::ranges::range_value_t<TSourceCont>, TSize >);

        calculate_node_weights_bf(graph, ignored_edges, ignored_nodes, check_abort);
        if (_neg_cycle) {
            return;
        }
        calculate_edge_weights(graph);

        auto reweighted_graph = copy_graph(graph, _edge_weights);
        if (threads == 1) {
            pathfind_multi(reweighted_graph, sources, ignored_edges, ignored_nodes, check_abort);
        } else {
            pathfind_multi_omp(reweighted_graph, sources, ignored_edges, ignored_nodes, check_abort, threads);
        }
    }

    /// @brief Attempts to compute the all pairs shortest paths for the given graph.
    ///
    /// If a negative cycle is detected, an edge will be removed from the graph based on the
    /// given edge scores and culling order.
    ///
    /// @param graph The graph to search.
    /// @param edge_scores The scores for each edge to be used when removing edges in negative cycles.
    /// @param cull_ascending Whether culling is performed in ascending or descending order of score.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @param threads The maximum number of threads to use during pathfinding. If <= 0, will be set to
    /// the number of processors available.
    void search_all_remove_cycles(
        const GraphViewType& graph,
        const WeightVec& edge_scores, bool cull_ascending,
        const CheckAbortFunc& check_abort = CheckAbortNoop,
        int threads = 0
    ) {
        search_all_remove_cycles(graph, {}, {}, edge_scores, cull_ascending, check_abort, threads);
    }

    /// @brief Attempts to compute the all pairs shortest paths for the given graph.
    ///
    /// If a negative cycle is detected, an edge will be removed from the graph based on the
    /// given edge scores and culling order.
    ///
    /// @param graph The graph to search.
    /// @param ignored_edges Set of edge IDs to ignore during pathfinding.
    /// @param ignored_nodes Set of node IDs to ignore during pathfinding.
    /// @param edge_scores The scores for each edge to be used when removing edges in negative cycles.
    /// @param cull_ascending Whether culling is performed in ascending or descending order of score.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @param threads The maximum number of threads to use during pathfinding. If <= 0, will be set to
    /// the number of processors available.
    void search_all_remove_cycles(
        const GraphViewType& graph,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const WeightVec& edge_scores, bool cull_ascending,
        const CheckAbortFunc& check_abort = CheckAbortNoop,
        int threads = 0
    ) {
        calculate_node_weights_ibf(graph, edge_scores, cull_ascending, ignored_edges, ignored_nodes, check_abort);
        calculate_edge_weights(graph);

        // Add any removed edges to the set of ignored edges during pathfinding
        EdgeIdSet all_ignored_edges(ignored_edges);
        for (auto edge_id : _removed_edges) {
            all_ignored_edges.emplace(edge_id);
        }

        auto reweighted_graph = copy_graph(graph, _edge_weights);
        if (threads == 1) {
            pathfind_all(reweighted_graph, all_ignored_edges, ignored_nodes, check_abort);
        } else {
            pathfind_all_omp(reweighted_graph, all_ignored_edges, ignored_nodes, check_abort, threads);
        }
    }

    /// @brief Attempts to compute the shortest paths for the given graph from the specified
    /// source nodes.
    ///
    /// If a negative cycle is detected, an edge will be removed from the graph based on the
    /// given edge scores and culling order.
    ///
    /// @param graph The graph to search.
    /// @param sources A container of source node IDs to pathfind from.
    /// @param edge_scores The scores for each edge to be used when removing edges in negative cycles.
    /// @param cull_ascending Whether culling is performed in ascending or descending order of score.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @param threads The maximum number of threads to use during pathfinding. If <= 0, will be set to
    /// the number of processors available.
    /// @tparam TSourceCont A container of source IDs.
    template <typename TSourceCont>
    void search_some_remove_cycles(
        const GraphViewType& graph, const TSourceCont& sources,
        const WeightVec& edge_scores, bool cull_ascending,
        const CheckAbortFunc& check_abort = CheckAbortNoop,
        int threads = 0
    ) {
        search_some_remove_cycles(graph, sources, {}, {}, edge_scores, cull_ascending, check_abort, threads);
    }

    /// @brief Attempts to compute the shortest paths for the given graph from the specified
    /// source nodes.
    ///
    /// If a negative cycle is detected, an edge will be removed from the graph based on the
    /// given edge scores and culling order.
    ///
    /// @param graph The graph to search.
    /// @param sources A container of source node IDs to pathfind from.
    /// @param ignored_edges Set of edge IDs to ignore during pathfinding.
    /// @param ignored_nodes Set of node IDs to ignore during pathfinding.
    /// @param edge_scores The scores for each edge to be used when removing edges in negative cycles.
    /// @param cull_ascending Whether culling is performed in ascending or descending order of score.
    /// @param check_abort Function that should throw an exception if execution should be aborted.
    /// @param threads The maximum number of threads to use during pathfinding. If <= 0, will be set to
    /// the number of processors available.
    /// @tparam TSourceCont A container of source IDs.
    template <typename TSourceCont>
    void search_some_remove_cycles(
        const GraphViewType& graph, const TSourceCont& sources,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const WeightVec& edge_scores, bool cull_ascending,
        const CheckAbortFunc& check_abort = CheckAbortNoop,
        int threads = 0
    ) {
        static_assert(std::ranges::range<TSourceCont>);
        static_assert(std::convertible_to< std::ranges::range_value_t<TSourceCont>, TSize >);

        calculate_node_weights_ibf(graph, edge_scores, cull_ascending, ignored_edges, ignored_nodes, check_abort);
        calculate_edge_weights(graph);

        // Add any removed edges to the set of ignored edges during pathfinding
        EdgeIdSet all_ignored_edges(ignored_edges);
        for (auto edge_id : _removed_edges) {
            all_ignored_edges.emplace(edge_id);
        }

        auto reweighted_graph = copy_graph(graph, _edge_weights);
        if (threads == 1) {
            pathfind_multi(reweighted_graph, sources, all_ignored_edges, ignored_nodes, check_abort);
        } else {
            pathfind_multi_omp(reweighted_graph, sources, all_ignored_edges, ignored_nodes, check_abort, threads);
        }
    }

private:
    GraphType copy_graph(const GraphViewType& graph) {
        bool weighted = graph.IsWeighted();
        GraphType copy;
        for (const auto& node : graph.Nodes()) {
            copy.CreateNode(node.id);
        }
        for (const auto& edge : graph.Edges()) {
            double weight = weighted ? graph.GetWeight(edge.id) : 1.0;
            copy.CreateEdge(
                edge.from, edge.to,
                mg_graph::GraphType::kDirectedGraph, edge.id, true, weight
            );
        }
        return copy;
    }

    GraphType copy_graph(const GraphViewType& graph, const WeightVec& weights) {
        GraphType copy;
        for (const auto& node : graph.Nodes()) {
            copy.CreateNode(node.id);
        }
        for (const auto& edge : graph.Edges()) {
            double weight = weights[edge.id];
            copy.CreateEdge(
                edge.from, edge.to,
                mg_graph::GraphType::kDirectedGraph, edge.id, true, weight
            );
        }
        return copy;
    }

    void calculate_node_weights_bf(
        const GraphViewType& graph,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort
    ) {
        _node_weights.clear();
        _neg_cycle.reset();
        _removed_edges.clear();

        // Create a copy of the graph with a new node that has 0-weight directed edges to all other nodes.
        auto copy = copy_graph(graph);
        TSize phantom_node_id = copy.Nodes().size();

        copy.CreateNode(phantom_node_id);
        for (TSize node_id = 0; node_id < phantom_node_id; node_id++) {
            copy.CreateEdge(phantom_node_id, node_id, mg_graph::GraphType::kDirectedGraph, std::nullopt, true, 0.0);
        }

        BellmanFordPF pathfinder;
        pathfinder.search(copy, phantom_node_id, ignored_edges, ignored_nodes, check_abort);

        if (pathfinder.has_negative_cycle()) {
            _neg_cycle = pathfinder.negative_cycle();
            return;
        }

        const auto& distances = pathfinder.distances();
        _node_weights.assign(phantom_node_id, 0.0);
        for (size_t i = 0; i < phantom_node_id; i++) {
            _node_weights[i] = distances[i];
        }
    }

    void calculate_node_weights_ibf(
        const GraphViewType& graph,
        const WeightVec& scores, bool cull_ascending,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort
    ) {
        _node_weights.clear();
        _neg_cycle.reset();
        _removed_edges.clear();

        // Create a copy of the graph with a new node that has 0-weight directed edges to all other nodes.
        auto copy = copy_graph(graph);
        TSize phantom_node_id = copy.Nodes().size();
        // Need to copy the scores and add entries for our new edges
        WeightVec all_scores;
        all_scores.reserve(scores.size() + phantom_node_id);
        for (auto score : scores) {
            all_scores.push_back(score);
        }

        copy.CreateNode(phantom_node_id);
        for (TSize node_id = 0; node_id < phantom_node_id; node_id++) {
            copy.CreateEdge(phantom_node_id, node_id, mg_graph::GraphType::kDirectedGraph, std::nullopt, true, 0.0);
            all_scores.push_back(0.0); // Doesn't matter what we add here really
        }

        IteraveBellmanFordPF pathfinder;
        pathfinder.remove_cycles(copy, phantom_node_id, ignored_edges, ignored_nodes, all_scores, cull_ascending, check_abort);

        for (const auto& edge_id : pathfinder.removed_edges()) {
            _removed_edges.emplace(edge_id);
        }

        const auto& distances = pathfinder.distances();
        _node_weights.assign(phantom_node_id, 0.0);
        for (size_t i = 0; i < phantom_node_id; i++) {
            _node_weights[i] = distances[i];
        }
    }

    void calculate_edge_weights(
        const GraphViewType& graph
    ) {
        const auto& edges = graph.Edges();

        WeightVec orig_edge_weights(edges.size(), 1.0);
        if (graph.IsWeighted()) {
            for (const auto& edge : graph.Edges()) {
                orig_edge_weights[edge.id] = graph.GetWeight(edge.id);
            }
        }

        _edge_weights.assign(edges.size(), 0.0);
        for (const auto& edge : graph.Edges()) {
            _edge_weights[edge.id] = orig_edge_weights[edge.id] + _node_weights[edge.from] - _node_weights[edge.to];
        }
    }

    void pathfind_all(
        const GraphViewType& graph,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort
    ) {
        size_t num_nodes = graph.Nodes().size();
        _pathfinders.clear();
        for (TSize source_id = 0; source_id < num_nodes; source_id++) {
            _pathfinders.emplace_back(std::move(std::make_unique<DijkstraPF>()));
            _pathfinders[source_id]->search(graph, source_id, ignored_edges, ignored_nodes, check_abort);
        }
    }

    void pathfind_all_omp(
        const GraphViewType& graph,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort, int threads
    ) {
        if (threads <=0) {
            threads = std::thread::hardware_concurrency();
        }

        size_t num_nodes = graph.Nodes().size();
        _pathfinders.clear();
        for (size_t source_id = 0; source_id < num_nodes; source_id++) {
            _pathfinders.emplace_back(std::move(std::make_unique<DijkstraPF>()));
        }

        omp_set_dynamic(threads);
        #pragma omp parallel for
        for (TSize source_id = 0; source_id < num_nodes; source_id++) {
            _pathfinders[source_id]->search(graph, source_id, ignored_edges, ignored_nodes, check_abort);
        }
    }

    template <typename TSourceCont>
    void pathfind_multi(
        const GraphViewType& graph, const TSourceCont& sources,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort
    ) {
        static_assert(std::ranges::range<TSourceCont>);
        static_assert(std::convertible_to< std::ranges::range_value_t<TSourceCont>, TSize >);

        size_t num_nodes = graph.Nodes().size();
        for (const auto& source_id : sources) {
            if (source_id >= num_nodes) {
                throw std::invalid_argument(
                    fmt::format("Invalid source node ID {}, should be in range [0, {}]", source_id, num_nodes)
                );
            }
        }

        _pathfinders.clear();
        for (size_t source_id = 0; source_id < num_nodes; source_id++) {
            _pathfinders.emplace_back(nullptr);
        }

        for (const auto& source_id : sources) {
            // Check if already initialized, and thus pathfound.
            if (_pathfinders[source_id]) continue;

            _pathfinders[source_id] = std::move(std::make_unique<DijkstraPF>());
            _pathfinders[source_id]->search(graph, source_id, ignored_edges, ignored_nodes, check_abort);
        }
    }

    template <typename TSourceCont>
    void pathfind_multi_omp(
        const GraphViewType& graph, const TSourceCont& sources,
        const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
        const CheckAbortFunc& check_abort, int threads
    ) {
        static_assert(std::ranges::range<TSourceCont>);
        static_assert(std::convertible_to< std::ranges::range_value_t<TSourceCont>, TSize >);

        if (threads <=0) {
            threads = std::thread::hardware_concurrency();
        }

        size_t num_nodes = graph.Nodes().size();
        for (const auto& source_id : sources) {
            if (source_id >= num_nodes) {
                throw std::invalid_argument(
                    fmt::format("Invalid source node ID {}, should be in range [0, {}]", source_id, num_nodes)
                );
            }
        }

        _pathfinders.clear();
        for (size_t source_id = 0; source_id < num_nodes; source_id++) {
            _pathfinders.emplace_back(nullptr);
        }

        NodeIdVec sources_to_search;
        for (const auto& source_id : sources) {
            // Check if already initialized
            if (_pathfinders[source_id]) continue;

            _pathfinders[source_id] = std::move(std::make_unique<DijkstraPF>());
            sources_to_search.push_back(source_id);
        }

        omp_set_dynamic(threads);
        #pragma omp parallel for
        for (auto source_id : sources_to_search) {
            _pathfinders[source_id]->search(graph, source_id, ignored_edges, ignored_nodes, check_abort);
        }
    }
};

} // namespace shortest_paths