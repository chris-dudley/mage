#include <mgp.hpp>
#include <mg_utils.hpp>

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <ranges>
#include <optional>

#include "algorithm/shortest_path.hpp"
#include "algorithm/dijkstra.hpp"
#include "algorithm/yens.hpp"
#include "algorithm/procedures.hpp"
#include "algorithm/bellman_ford.hpp"
#include "algorithm/iterative_bf.hpp"
#include "algorithm/johnsons.hpp"
#include "algorithm/disjoint.hpp"

namespace {

constexpr const int64_t ALGO_DIJKSTRA = 0;

// For whatever reason, calling ->GetInnerNodeId directly on the unique_ptr returned by Get(Sub)graphView breaks,
// and will throw an mg_exception::InvalidIDException if the internal memgraph node ID is outside the range of
// the internal IDs in the view. (e.g. there are 8 nodes in the view and you're trying to lookup memgraph_id 9)
// Creating a function that takes a const mg_graph::GraphView<>& and passing in *ptr makes it work properly, and
// I'm not sure why.
uint64_t GetInnerNodeId(const mg_graph::GraphView<>& graph, uint64_t memgraph_id) {
    return graph.GetInnerNodeId(memgraph_id);
}

std::optional<uint64_t> GetInnerNodeIdOpt(const mg_graph::GraphView<>& graph, uint64_t memgraph_id) {
    if (graph.NodeExists(memgraph_id)) {
        return graph.GetInnerNodeId(memgraph_id);
    }
    return std::nullopt;
}

uint64_t GetInnerEdgeId(const mg_graph::GraphView<>& graph, uint64_t memgraph_id) {
    return graph.GetInnerEdgeId(memgraph_id);
}

// Add wrappers for GetMemgraph(Node|Edge)Id just in case.
uint64_t GetMemgraphNodeId(const mg_graph::GraphView<>& graph, uint64_t inner_id) {
    return graph.GetMemgraphNodeId(inner_id);
}
uint64_t GetMemgraphEdgeId(const mg_graph::GraphView<>& graph, uint64_t inner_id) {
    return graph.GetMemgraphEdgeId(inner_id);
}

mgp::Path TranslatePath(const mgp::Graph& graph, const mg_graph::GraphView<>& graph_view, const shortest_paths::Path<>& path) {
    auto source_mgid = mgp::Id::FromUint(GetMemgraphNodeId(graph_view, path.nodes[0]));
    mgp::Node source_node = graph.GetNodeById(source_mgid);
    mgp::Path result_path(source_node);

    mgp::Node current_node = source_node;
    for (size_t edge_index = 0; edge_index < path.size(); edge_index++) {
        auto from_node_mgid = GetMemgraphNodeId(graph_view, path.nodes[edge_index]);
        auto to_node_mgid = GetMemgraphNodeId(graph_view, path.nodes[edge_index+1]);
        auto edge_mgid = GetMemgraphEdgeId(graph_view, path.edges[edge_index]);

        if (from_node_mgid != current_node.Id().AsUint()) {
            throw std::logic_error("From node in edge does not match current node!");
        }

        bool found_edge = false;
        // Search for the matching edge since there's no easy way of pulling out a specific
        // edge from the graph by ID.
        for (const auto& relationship : current_node.OutRelationships()) {
            if (relationship.Id().AsUint() != edge_mgid) {
                continue;
            }
            auto next_node = relationship.To();
            if (next_node.Id().AsUint() != to_node_mgid) {
                throw std::logic_error("To node in edge does not match to node in relationship!");
            }

            found_edge = true;

            result_path.Expand(relationship);
            current_node = next_node;
            break;
        }

        if (!found_edge) {
            throw std::logic_error("Unable to find matching out relationship");
        }
    }

    return result_path;
}

mgp::List TranslateNodeList(const mgp::Graph& graph, const mg_graph::GraphView<>& graph_view, const std::vector<uint64_t> node_ids) {
    mgp::List result(node_ids.size());
    for (auto node_id : node_ids) {
        auto mgid = GetMemgraphNodeId(graph_view, node_id);
        result.Append(mgp::Value(graph.GetNodeById(mgp::Id::FromUint(mgid))));
    }
    return result;
}

mgp::List TranslateEdgeList(const mgp::Graph& graph, const mg_graph::GraphView<>& graph_view, const std::vector<uint64_t>& edge_ids) {
    std::vector<mgp::Value> edges(edge_ids.size());
    std::unordered_multimap<uint64_t, size_t> mgid_to_index;
    for (size_t index = 0; index < edge_ids.size(); index++) {
        auto mgid = GetMemgraphEdgeId(graph_view, edge_ids[index]);
        mgid_to_index.emplace(mgid, index);
    }
    for (auto edge : graph.Relationships()) {
        auto range = mgid_to_index.equal_range(edge.Id().AsUint());
        for (auto it = range.first; it != range.second; it++) {
            edges[it->second] = mgp::Value(edge);
        }
    }

    return mgp::List(edges);
}

template<typename TRange>
mgp::List ToList(TRange&& range) {
    mgp::List result(std::ranges::size(range));
    for (auto&& value : range) {
        result.Append(mgp::Value(value));
    }
    return result;
}

void YensKShortestPaths(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        Source, Sink, K, WeightProp, DefaultWeight, Threads
    };

    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        const auto source_node = arguments[ArgIdx::Source].ValueNode();
        const auto sink_node = arguments[ArgIdx::Sink].ValueNode();
        const auto K = arguments[ArgIdx::K].ValueInt();
        const auto weight_property = arguments[ArgIdx::WeightProp].ValueString();
        const auto default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();
        const int threads = arguments[ArgIdx::Threads].ValueInt();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.data() : nullptr;

        if (K < 0) {
            throw mgp::ValueException("K cannot be negative");
        }
        if (K == 0) {
            // Why did you even call this procedure?
            return;
        }

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );

        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        const auto source_mgid = source_node.Id().AsUint();
        const auto sink_mgid = sink_node.Id().AsUint();

        const auto maybe_source_id = GetInnerNodeIdOpt(graph_view, source_mgid);
        const auto maybe_sink_id = GetInnerNodeIdOpt(graph_view, sink_mgid);
        if (!maybe_source_id || !maybe_sink_id) {
            // Either source or sink is not in subgraph, so there can be no paths between them.
            return;
        }
        const auto source_id = maybe_source_id.value();
        const auto sink_id = maybe_sink_id.value();

        auto abort_func = [&graph] () { graph.CheckMustAbort(); };
        const auto paths = shortest_paths::KShortestPaths(graph_view, source_id, sink_id, K, abort_func, threads);

        for (std::size_t path_index = 0; path_index < paths.size(); path_index++) {
            const auto& path = paths[path_index];

            mgp::Path result_path = TranslatePath(graph, graph_view, path);
            mgp::Node current_node = source_node;
            mgp::List costs(path.costs.size()); // accumulated cost at each node
            for (auto cost : path.costs) {
                costs.Append(mgp::Value(cost));
            }

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnIndex).c_str(), int64_t(path_index));
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), source_node);
            record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), sink_node);
            record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), path.total_cost);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
        }
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}


void BellmanFordProcedure(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        Source, Targets, WeightProp, DefaultWeight
    };

    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {

        const auto source_node = arguments[ArgIdx::Source].ValueNode();
        const auto target_nodes = arguments[ArgIdx::Targets].ValueList();
        const auto weight_property = arguments[ArgIdx::WeightProp].ValueString();
        const auto default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.data() : nullptr;

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );

        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        const auto source_mgid = source_node.Id().AsUint();
        const auto maybe_source_id = GetInnerNodeIdOpt(graph_view, source_mgid);
        if (!maybe_source_id) {
            // Source is not in subgraph, so there can be no paths from it.
            return;
        }
        const auto source_id = maybe_source_id.value();

        shortest_paths::BellmanFordPathfinder<uint64_t> pathfinder(graph_view, source_id);
        if (pathfinder.has_negative_cycle()) {
            auto cycle = pathfinder.negative_cycle().value();
            mgp::Path result_path = TranslatePath(graph, graph_view, cycle);
            mgp::Node cycle_source = result_path.GetNodeAt(0);
            mgp::List costs(cycle.costs.size()); // accumulated cost at each node
            for (auto cost : cycle.costs) {
                costs.Append(mgp::Value(cost));
            }

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), true);
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), cycle_source);
            record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), cycle_source);
            record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), cycle.total_cost);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
            return;
        }

        std::vector<uint64_t> targets;
        if (target_nodes.Empty()) {
            // If no targets specified, return paths to all reachable nodes
            targets = pathfinder.reachable_nodes();
        } else {
            for (auto target_nodes_elem : target_nodes ) {
                auto target_node = target_nodes_elem.ValueNode();
                auto node_id = GetInnerNodeId(graph_view, target_node.Id().AsUint());
                if (pathfinder.has_path_to(node_id)) {
                    targets.push_back(node_id);
                }
            }
        }

        for (auto target_id : targets) {
            auto target_mgid = GetMemgraphNodeId(graph_view, target_id);

            auto path = pathfinder.path_to(target_id);
            if (path.empty()) continue;

            mgp::Node target_node = graph.GetNodeById(mgp::Id::FromUint(target_mgid));
            mgp::Path result_path = TranslatePath(graph, graph_view, path);
            mgp::List costs(path.costs.size()); // accumulated cost at each node
            for (auto cost : path.costs) {
                costs.Append(mgp::Value(cost));
            }

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), false);
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), source_node);
            record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), target_node);
            record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), path.total_cost);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
        }
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}

std::vector<double> GetScores(const mgp::Graph& graph, const mg_graph::GraphView<>& graph_view, const std::string& score_prop_name, double default_score) {
    const auto& edges = graph_view.Edges();
    size_t num_edges = edges.size();
    std::vector<double> scores(num_edges, default_score);
    if (score_prop_name.empty() && !graph_view.IsWeighted()) {
        return scores;
    }
    if (score_prop_name.empty()) {
        // Default to edge weights
        for (uint64_t edge_id = 0; edge_id < num_edges; edge_id++) {
            scores[edge_id] = graph_view.GetWeight(edge_id);
        }
        return scores;
    }
    for (const auto& relationship : graph.Relationships()) {
        auto inner_id = GetInnerEdgeId(graph_view, relationship.Id().AsUint());
        auto score_prop = relationship.GetProperty(score_prop_name);
    
        // To match the handling of weights in mg_graph, properties that aren't doubles or ints will receive
        // the default score.
        if (score_prop.IsNumeric()) {
            scores[inner_id] = score_prop.ValueNumeric();
        }
    }

    return scores;
}

void IterativeBellmanFordProcedure(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        Source, Target, WeightProp, ScoreProp, DefaultWeight, DefaultScore, CullAscending
    };

    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        const auto source_node = arguments[ArgIdx::Source].ValueNode();
        const auto target_node = arguments[ArgIdx::Target].ValueNode();
        const auto weight_property = std::string(arguments[ArgIdx::WeightProp].ValueString());
        const auto score_property = std::string(arguments[ArgIdx::ScoreProp].ValueString());
        double default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();
        double default_score = arguments[ArgIdx::DefaultScore].ValueDouble();
        bool cull_ascending = arguments[ArgIdx::CullAscending].ValueBool();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.c_str() : nullptr;

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );
        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        std::vector<double> scores = GetScores(graph, graph_view, score_property, default_score);


        const auto source_mgid = source_node.Id().AsUint();
        const auto target_mgid = target_node.Id().AsUint();
        const auto maybe_source_id = GetInnerNodeIdOpt(graph_view, source_mgid);
        const auto maybe_target_id = GetInnerNodeIdOpt(graph_view, target_mgid);
        if (!maybe_source_id || !maybe_target_id) {
            // Either source or target is not in subgraph, so there can be no paths between them.
            return;
        }
        const auto source_id = maybe_source_id.value();
        const auto target_id = maybe_target_id.value();
        auto abort_func = [&graph] () { graph.CheckMustAbort(); };

        shortest_paths::IterativeBellmanFordPathfinder<uint64_t> pathfinder(scores, cull_ascending);
        auto path = pathfinder.search(graph_view, source_id, target_id, {}, {}, abort_func);
        if (path.empty()) {
            return;
        }

        mgp::Path result_path = TranslatePath(graph, graph_view, path);
        mgp::List costs(path.costs.size()); // accumulated cost at each node
        for (auto cost : path.costs) {
            costs.Append(mgp::Value(cost));
        }

        auto record = record_factory.NewRecord();
        record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), source_node);
        record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), target_node);
        record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), path.total_cost);
        record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
        record.Insert(
            std::string(shortest_paths::kReturnNumEdgesRemoved).c_str(),
            static_cast<int64_t>(pathfinder.num_edges_removed())
        );
        record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}

void IBF_SSSP(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        Source, Targets, WeightProp, ScoreProp, DefaultWeight, DefaultScore, CullAscending
    };
    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        const auto source_node = arguments[ArgIdx::Source].ValueNode();
        const auto target_nodes = arguments[ArgIdx::Targets].ValueList();
        const auto weight_property = std::string(arguments[ArgIdx::WeightProp].ValueString());
        const auto score_property = std::string(arguments[ArgIdx::ScoreProp].ValueString());
        double default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();
        double default_score = arguments[ArgIdx::DefaultScore].ValueDouble();
        bool cull_ascending = arguments[ArgIdx::CullAscending].ValueBool();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.c_str() : nullptr;

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );
        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        std::vector<double> scores = GetScores(graph, graph_view, score_property, default_score);

        const auto source_mgid = source_node.Id().AsUint();
        const auto maybe_source_id = GetInnerNodeIdOpt(graph_view, source_mgid);
        if (!maybe_source_id) {
            // Source is not in subgraph, so there can be no paths from it.
            return;
        }
        const auto source_id = maybe_source_id.value();
        auto abort_func = [&graph] () { graph.CheckMustAbort(); };

        shortest_paths::IterativeBellmanFordPathfinder<uint64_t> pathfinder(scores, cull_ascending);
        pathfinder.remove_cycles(graph_view, source_id, {}, {}, abort_func);

        std::vector<uint64_t> targets;
        if (target_nodes.Empty()) {
            targets = pathfinder.reachable_nodes();
        } else {
            for (auto target_nodes_elem : target_nodes) {
                auto target_node = target_nodes_elem.ValueNode();
                auto node_id = GetInnerNodeId(graph_view, target_node.Id().AsUint());
                if (pathfinder.has_path_to(node_id)) {
                    targets.push_back(node_id);
                }
            }
        }

        for (auto target_id : targets) {
            if (target_id == source_id) continue;

            auto target_mgid = GetMemgraphNodeId(graph_view, target_id);

            auto path = pathfinder.path_to(target_id);
            if (path.empty()) continue;

            mgp::Node target_node = graph.GetNodeById(mgp::Id::FromUint(target_mgid));
            mgp::Path result_path = TranslatePath(graph, graph_view, path);
            mgp::List costs(path.costs.size()); // accumulated cost at each node
            for (auto cost : path.costs) {
                costs.Append(mgp::Value(cost));
            }

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), source_node);
            record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), target_node);
            record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), path.total_cost);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
        }
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}


void IBF_SSSP_Subgraph(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        Source, WeightProp, ScoreProp, DefaultWeight, DefaultScore, CullAscending
    };
    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        const auto source_node = arguments[ArgIdx::Source].ValueNode();
        const auto weight_property = std::string(arguments[ArgIdx::WeightProp].ValueString());
        const auto score_property = std::string(arguments[ArgIdx::ScoreProp].ValueString());
        double default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();
        double default_score = arguments[ArgIdx::DefaultScore].ValueDouble();
        bool cull_ascending = arguments[ArgIdx::CullAscending].ValueBool();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.c_str() : nullptr;

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );
        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        std::vector<double> scores = GetScores(graph, graph_view, score_property, default_score);

        const auto source_mgid = source_node.Id().AsUint();
        const auto maybe_source_id = GetInnerNodeIdOpt(graph_view, source_mgid);
        if (!maybe_source_id) {
            // Source not in (sub)graph, so resulting subgraph is empty.
            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNodes).c_str(), mgp::List());
            record.Insert(std::string(shortest_paths::kReturnEdges).c_str(), mgp::List());
            record.Insert(std::string(shortest_paths::kReturnNumEdgesRemoved).c_str(), 0L);
            record.Insert(std::string(shortest_paths::kReturnRemovedEdges).c_str(), mgp::List());
            return;
        }

        const auto source_id = maybe_source_id.value();
        auto abort_func = [&graph] () { graph.CheckMustAbort(); };

        shortest_paths::IterativeBellmanFordPathfinder<uint64_t> pathfinder(scores, cull_ascending);
        pathfinder.remove_cycles(graph_view, source_id, {}, {}, abort_func);

        // Add all reachable nodes to result set
        const auto reachable_node_ids = pathfinder.reachable_nodes();
        mgp::List result_nodes(reachable_node_ids.size());
        for (auto node_id : reachable_node_ids) {
            auto mgid = mgp::Id::FromUint(GetMemgraphNodeId(graph_view, node_id));
            result_nodes.Append(mgp::Value(graph.GetNodeById(mgid)));
        }

        // Set up sets of memgraph IDs for included and removed edges, so we can more easily check for inclusion
        // in either set while iterating over the edges.
        auto edge_ids_used = pathfinder.edges_used();
        auto removed_edge_ids = pathfinder.removed_edges();
        std::unordered_set<uint64_t> edges_used_mgid_set;
        std::unordered_set<uint64_t> removed_edges_mgid_set;
        edges_used_mgid_set.reserve(edge_ids_used.size());
        removed_edges_mgid_set.reserve(removed_edge_ids.size());

        for (auto edge_id : edge_ids_used) {
            edges_used_mgid_set.emplace(GetMemgraphEdgeId(graph_view, edge_id));
        }
        for (auto edge_id : removed_edge_ids) {
            removed_edges_mgid_set.emplace(GetMemgraphEdgeId(graph_view, edge_id));
        }

        mgp::List result_edges(edges_used_mgid_set.size());
        mgp::List result_removed_edges(removed_edges_mgid_set.size());

        for (const auto& edge : graph.Relationships()) {
            auto edge_mgid = edge.Id().AsUint();
            if (edges_used_mgid_set.contains(edge_mgid)) {
                result_edges.Append(mgp::Value(edge));
            } else if (removed_edges_mgid_set.contains(edge_mgid)) {
                result_removed_edges.Append(mgp::Value(edge));
            }
        }

        auto record = record_factory.NewRecord();
        record.Insert(std::string(shortest_paths::kReturnNodes).c_str(), result_nodes);
        record.Insert(std::string(shortest_paths::kReturnEdges).c_str(), result_edges);
        record.Insert(std::string(shortest_paths::kReturnNumEdgesRemoved).c_str(), static_cast<int64_t>(removed_edges_mgid_set.size()));
        record.Insert(std::string(shortest_paths::kReturnRemovedEdges).c_str(), result_removed_edges);
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}

void Johnsons_Subgraph(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        WeightProp, ScoreProp, DefaultWeight, DefaultScore, CullAscending, Threads
    };
    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        const auto weight_property = std::string(arguments[ArgIdx::WeightProp].ValueString());
        const auto score_property = std::string(arguments[ArgIdx::ScoreProp].ValueString());
        double default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();
        double default_score = arguments[ArgIdx::DefaultScore].ValueDouble();
        bool cull_ascending = arguments[ArgIdx::CullAscending].ValueBool();
        int threads = arguments[ArgIdx::Threads].ValueInt();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.c_str() : nullptr;

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );
        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        auto abort_func = [&graph] () { graph.CheckMustAbort(); };

        shortest_paths::JohnsonsPathfinder<> pathfinder;
        if (!score_property.empty()) {
            std::vector<double> scores = GetScores(graph, graph_view, score_property, default_score);
            pathfinder.search_all_remove_cycles(graph_view, scores, cull_ascending, abort_func, threads);
        } else {
            pathfinder.search_all(graph_view, abort_func, threads);
        }

        if (pathfinder.has_negative_cycle()) {
            auto cycle = pathfinder.negative_cycle().value();
            mgp::List nodes(cycle.nodes.size());
            mgp::List edges(cycle.edges.size());
            // Empty list return parameters
            mgp::List removed_edges;
            mgp::List edge_weights;

            mgp::Path result_path = TranslatePath(graph, graph_view, cycle);
            mgp::Node cycle_source = result_path.GetNodeAt(0);
            mgp::List costs(cycle.costs.size()); // accumulated cost at each node
            for (auto cost : cycle.costs) {
                costs.Append(mgp::Value(cost));
            }
            nodes.Append(mgp::Value(cycle_source));
            for (size_t edge_index = 0; edge_index < cycle.size(); edge_index++) {
                edges.Append(mgp::Value(result_path.GetRelationshipAt(edge_index)));
                nodes.Append(mgp::Value(result_path.GetNodeAt(edge_index+1)));
            }

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), true);
            record.Insert(std::string(shortest_paths::kReturnNodes).c_str(), nodes);
            record.Insert(std::string(shortest_paths::kReturnEdges).c_str(), edges);
            record.Insert(std::string(shortest_paths::kReturnNodeWeights).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnEdgeWeights).c_str(), edge_weights);
            record.Insert(std::string(shortest_paths::kReturnNumEdgesRemoved).c_str(), 0L);
            record.Insert(std::string(shortest_paths::kReturnRemovedEdges).c_str(), removed_edges);
            return;
        }

        // List of nodes in same order as view, so the edge weights match.
        mgp::List nodes(graph_view.Nodes().size());
        for (const auto& node : graph_view.Nodes()) {
            auto node_id = mgp::Id::FromUint(GetMemgraphNodeId(graph_view, node.id));
            nodes.Append(mgp::Value(graph.GetNodeById(node_id)));
        }
        auto node_weights = ToList(pathfinder.node_weights());
        const auto& inner_edge_weights = pathfinder.edge_weights();
        std::unordered_map<uint64_t, double> mgid_to_edge_weight;

        int64_t num_edges_removed = static_cast<int64_t>(pathfinder.num_edges_removed());

        std::vector<bool> checked_edges(graph_view.Edges().size(), false);
        std::unordered_set<uint64_t> edges_used_mgid_set;
        std::unordered_set<uint64_t> removed_edges_mgid_set;
        removed_edges_mgid_set.reserve(pathfinder.num_edges_removed());

        // Build set of edges used by any shortest path
        for (uint64_t source_id = 0; source_id < graph_view.Nodes().size(); source_id++) {
            for (auto edge_id : pathfinder.edges_used_from(source_id)) {
                if (checked_edges[edge_id]) {
                    continue;
                }
                checked_edges[edge_id] = true;
                auto mgid = GetMemgraphEdgeId(graph_view, edge_id);
                edges_used_mgid_set.emplace(mgid);
                mgid_to_edge_weight[mgid] = inner_edge_weights[edge_id];
            }
        }
        for (auto edge_id : pathfinder.removed_edges()) {
            removed_edges_mgid_set.emplace(GetMemgraphEdgeId(graph_view, edge_id));
        }

        mgp::List edges(edges_used_mgid_set.size());
        mgp::List edge_weights(edges_used_mgid_set.size());
        mgp::List removed_edges(removed_edges_mgid_set.size());
        for (const auto& edge : graph.Relationships()) {
            auto edge_mgid = edge.Id().AsUint();
            if (edges_used_mgid_set.contains(edge_mgid)) {
                edges.Append(mgp::Value(edge));
                edge_weights.Append(mgp::Value(mgid_to_edge_weight[edge_mgid]));
            } else if (removed_edges_mgid_set.contains(edge_mgid)) {
                removed_edges.Append(mgp::Value(edge));
            }
        }

        auto record = record_factory.NewRecord();
        record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), false);
        record.Insert(std::string(shortest_paths::kReturnNodes).c_str(), nodes);
        record.Insert(std::string(shortest_paths::kReturnEdges).c_str(), edges);
        record.Insert(std::string(shortest_paths::kReturnNodeWeights).c_str(), node_weights);
        record.Insert(std::string(shortest_paths::kReturnEdgeWeights).c_str(), edge_weights);
        record.Insert(std::string(shortest_paths::kReturnNumEdgesRemoved).c_str(), num_edges_removed);
        record.Insert(std::string(shortest_paths::kReturnRemovedEdges).c_str(), removed_edges);
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}

void Johnsons_Source_Subgraphs(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        SourceNodes, WeightProp, ScoreProp, DefaultWeight, DefaultScore, CullAscending, Threads
    };
    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        auto source_nodes = arguments[ArgIdx::SourceNodes].ValueList();
        const auto weight_property = std::string(arguments[ArgIdx::WeightProp].ValueString());
        const auto score_property = std::string(arguments[ArgIdx::ScoreProp].ValueString());
        double default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();
        double default_score = arguments[ArgIdx::DefaultScore].ValueDouble();
        bool cull_ascending = arguments[ArgIdx::CullAscending].ValueBool();
        int threads = arguments[ArgIdx::Threads].ValueInt();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.c_str() : nullptr;

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );
        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        auto abort_func = [&graph] () { graph.CheckMustAbort(); };

        std::unordered_set<uint64_t> source_ids;
        for (const auto& value : source_nodes) {
            auto node = value.ValueNode();
            auto source_id = GetInnerNodeIdOpt(graph_view, node.Id().AsUint());
            // Ignore sources not in subgraph
            if (source_id) {
                source_ids.emplace(source_id.value());
            }
        }
        // If no nodes specified, search all
        if (source_ids.empty()) {
            for (auto node : graph_view.Nodes()) {
                source_ids.emplace(node.id);
            }
        }

        shortest_paths::JohnsonsPathfinder<> pathfinder;
        if (!score_property.empty()) {
            std::vector<double> scores = GetScores(graph, graph_view, score_property, default_score);
            pathfinder.search_some_remove_cycles(graph_view, source_ids, scores, cull_ascending, abort_func, threads);
        } else {
            pathfinder.search_some(graph_view, source_ids, abort_func, threads);
        }

        if (pathfinder.has_negative_cycle()) {
            auto cycle = pathfinder.negative_cycle().value();
            mgp::List nodes(cycle.nodes.size());
            mgp::List edges(cycle.edges.size());

            mgp::Path result_path = TranslatePath(graph, graph_view, cycle);
            mgp::Node cycle_source = result_path.GetNodeAt(0);
            auto costs = ToList(cycle.costs);

            nodes.Append(mgp::Value(cycle_source));
            for (size_t edge_index = 0; edge_index < cycle.size(); edge_index++) {
                edges.Append(mgp::Value(result_path.GetRelationshipAt(edge_index)));
                nodes.Append(mgp::Value(result_path.GetNodeAt(edge_index+1)));
            }

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), true);
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), cycle_source);
            record.Insert(std::string(shortest_paths::kReturnNodes).c_str(), nodes);
            record.Insert(std::string(shortest_paths::kReturnEdges).c_str(), edges);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            return;
        }

        for (auto source_id : source_ids) {
            auto source_mgid = GetMemgraphNodeId(graph_view, source_id);
            auto source = graph.GetNodeById(mgp::Id::FromUint(source_mgid));
            auto nodes = TranslateNodeList(graph, graph_view, pathfinder.reachable_nodes_from(source_id));
            auto edges = TranslateEdgeList(graph, graph_view, pathfinder.edges_used_from(source_id));
            auto costs = ToList(pathfinder.distances_from(source_id));

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), false);
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), source);
            record.Insert(std::string(shortest_paths::kReturnNodes).c_str(), nodes);
            record.Insert(std::string(shortest_paths::kReturnEdges).c_str(), edges);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
        }
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}

void Johnsons_Paths(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        SourceNodes, TargetNodes, WeightProp, ScoreProp, DefaultWeight, DefaultScore, CullAscending, Threads
    };
    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        auto source_nodes = arguments[ArgIdx::SourceNodes].ValueList();
        auto target_nodes = arguments[ArgIdx::TargetNodes].ValueList();
        const auto weight_property = std::string(arguments[ArgIdx::WeightProp].ValueString());
        const auto score_property = std::string(arguments[ArgIdx::ScoreProp].ValueString());
        double default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();
        double default_score = arguments[ArgIdx::DefaultScore].ValueDouble();
        bool cull_ascending = arguments[ArgIdx::CullAscending].ValueBool();
        int threads = arguments[ArgIdx::Threads].ValueInt();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.c_str() : nullptr;

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );
        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        auto abort_func = [&graph] () { graph.CheckMustAbort(); };

        std::unordered_set<uint64_t> source_ids;
        for (const auto& value : source_nodes) {
            auto node = value.ValueNode();
            auto source_id = GetInnerNodeIdOpt(graph_view, node.Id().AsUint());
            // Ignore sources not in subgraph
            if (source_id) {
                source_ids.emplace(source_id.value());
            }
        }
        // If no nodes specified, search all
        if (source_ids.empty()) {
            for (auto node : graph_view.Nodes()) {
                source_ids.emplace(node.id);
            }
        }

        std::unordered_set<uint64_t> target_ids;
        for (const auto& value : target_nodes) {
            auto node = value.ValueNode();
            auto target_id = GetInnerNodeIdOpt(graph_view, node.Id().AsUint());
            // Ignore targets not in subgraph
            if (target_id) {
                target_ids.emplace(target_id.value());
            }
        }
        // If no nodes specified, target all
        if (target_ids.empty()) {
            for (auto node : graph_view.Nodes()) {
                target_ids.emplace(node.id);
            }
        }

        shortest_paths::JohnsonsPathfinder<> pathfinder;
        if (!score_property.empty()) {
            std::vector<double> scores = GetScores(graph, graph_view, score_property, default_score);
            pathfinder.search_some_remove_cycles(graph_view, source_ids, scores, cull_ascending, abort_func, threads);
        } else {
            pathfinder.search_some(graph_view, source_ids, abort_func, threads);
        }

        if (pathfinder.has_negative_cycle()) {
            auto cycle = pathfinder.negative_cycle().value();
            mgp::Path result_path = TranslatePath(graph, graph_view, cycle);
            mgp::Node cycle_source = result_path.GetNodeAt(0);
            auto costs = ToList(cycle.costs);

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), true);
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), cycle_source);
            record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), cycle_source);
            record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), cycle.total_cost);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
            return;
        }

        for (auto source_id : source_ids) {
            auto source_mgid = GetMemgraphNodeId(graph_view, source_id);
            auto source = graph.GetNodeById(mgp::Id::FromUint(source_mgid));
            for (auto target_id : target_ids) {
                if (target_id == source_id || !pathfinder.has_path(source_id, target_id)) {
                    continue;
                }

                auto target_mgid = GetMemgraphNodeId(graph_view, target_id);
                auto target = graph.GetNodeById(mgp::Id::FromUint(target_mgid));

                auto path = pathfinder.get_path(source_id, target_id);
                auto result_path = TranslatePath(graph, graph_view, path);
                auto costs = ToList(path.costs);

                auto record = record_factory.NewRecord();
                record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), false);
                record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), source);
                record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), target);
                record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), path.total_cost);
                record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
                record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
            }
        }
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}

void Johnsons_K_Shortest(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    using PathVec = shortest_paths::JohnsonsPathfinder<>::PathVec;

    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        SourceNode, TargetNode, K, WeightProp, ScoreProp, DefaultWeight, DefaultScore, CullAscending, Threads
    };
    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        const auto source_node = arguments[ArgIdx::SourceNode].ValueNode();
        const auto target_node = arguments[ArgIdx::TargetNode].ValueNode();
        const auto K = arguments[ArgIdx::K].ValueInt();
        const auto weight_property = std::string(arguments[ArgIdx::WeightProp].ValueString());
        const auto score_property = std::string(arguments[ArgIdx::ScoreProp].ValueString());
        double default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();
        double default_score = arguments[ArgIdx::DefaultScore].ValueDouble();
        bool cull_ascending = arguments[ArgIdx::CullAscending].ValueBool();
        int threads = arguments[ArgIdx::Threads].ValueInt();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.c_str() : nullptr;

        if (K < 0) {
            throw mgp::ValueException("K cannot be negative");
        }
        if (K == 0) {
            // Why did you even call this procedure?
            return;
        }

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );
        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        const auto source_mgid = source_node.Id().AsUint();
        const auto target_mgid = target_node.Id().AsUint();
        const auto maybe_source_id = GetInnerNodeIdOpt(graph_view, source_mgid);
        const auto maybe_target_id = GetInnerNodeIdOpt(graph_view, target_mgid);
        if (!maybe_source_id || !maybe_target_id) {
            // Either source or target is not in subgraph, so there can be no paths between them.
            return;
        }
        const auto source_id = maybe_source_id.value();
        const auto target_id = maybe_target_id.value();

        auto abort_func = [&graph] () { graph.CheckMustAbort(); };

        shortest_paths::JohnsonsPathfinder<> pathfinder;
        PathVec paths;
        if (!score_property.empty()) {
            std::vector<double> scores = GetScores(graph, graph_view, score_property, default_score);
            paths = pathfinder.k_shortest_paths_remove_cycles(graph_view, source_id, target_id, K, scores, cull_ascending, abort_func, threads);
        } else {
            paths = pathfinder.k_shortest_paths(graph_view, source_id, target_id, K, abort_func, threads);
        }

        if (pathfinder.has_negative_cycle()) {
            auto cycle = pathfinder.negative_cycle().value();
            mgp::Path result_path = TranslatePath(graph, graph_view, cycle);
            mgp::Node cycle_source = result_path.GetNodeAt(0);
            auto costs = ToList(cycle.costs);

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), true);
            record.Insert(std::string(shortest_paths::kReturnIndex).c_str(), 0L);
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), cycle_source);
            record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), cycle_source);
            record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), cycle.total_cost);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
            return;
        }

        for (std::size_t path_index = 0; path_index < paths.size(); path_index++) {
            const auto& path = paths[path_index];

            mgp::Path result_path = TranslatePath(graph, graph_view, path);
            mgp::Node current_node = source_node;
            auto costs = ToList(path.costs);

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), false);
            record.Insert(std::string(shortest_paths::kReturnIndex).c_str(), int64_t(path_index));
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), source_node);
            record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), target_node);
            record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), path.total_cost);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
        }
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}

void Johnsons_Disjoint_K_Shortest(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    using PathVec = shortest_paths::JohnsonsPathfinder<>::PathVec;

    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        SourceNode, TargetNode, K, WeightProp, ScoreProp, DefaultWeight, DefaultScore, CullAscending
    };
    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        const auto source_node = arguments[ArgIdx::SourceNode].ValueNode();
        const auto target_node = arguments[ArgIdx::TargetNode].ValueNode();
        const auto K = arguments[ArgIdx::K].ValueInt();
        const auto weight_property = std::string(arguments[ArgIdx::WeightProp].ValueString());
        const auto score_property = std::string(arguments[ArgIdx::ScoreProp].ValueString());
        double default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();
        double default_score = arguments[ArgIdx::DefaultScore].ValueDouble();
        bool cull_ascending = arguments[ArgIdx::CullAscending].ValueBool();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.c_str() : nullptr;

        if (K < 0) {
            throw mgp::ValueException("K cannot be negative");
        }

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );
        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        const auto source_mgid = source_node.Id().AsUint();
        const auto target_mgid = target_node.Id().AsUint();
        const auto maybe_source_id = GetInnerNodeIdOpt(graph_view, source_mgid);
        const auto maybe_target_id = GetInnerNodeIdOpt(graph_view, target_mgid);
        if (!maybe_source_id || !maybe_target_id) {
            // Either source or target is not in subgraph, so there can be no paths between them.
            return;
        }
        const auto source_id = maybe_source_id.value();
        const auto target_id = maybe_target_id.value();

        auto abort_func = [&graph] () { graph.CheckMustAbort(); };

        shortest_paths::JohnsonsPathfinder<> pathfinder;
        PathVec paths;
        if (!score_property.empty()) {
            std::vector<double> scores = GetScores(graph, graph_view, score_property, default_score);
            paths = pathfinder.disjoint_k_shortest_paths_remove_cycles(graph_view, source_id, target_id, K, scores, cull_ascending, abort_func);
        } else {
            paths = pathfinder.disjoint_k_shortest_paths(graph_view, source_id, target_id, K, abort_func);
        }

        if (pathfinder.has_negative_cycle()) {
            auto cycle = pathfinder.negative_cycle().value();
            mgp::Path result_path = TranslatePath(graph, graph_view, cycle);
            mgp::Node cycle_source = result_path.GetNodeAt(0);
            auto costs = ToList(cycle.costs);

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), true);
            record.Insert(std::string(shortest_paths::kReturnIndex).c_str(), 0L);
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), cycle_source);
            record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), cycle_source);
            record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), cycle.total_cost);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
            return;
        }

        for (std::size_t path_index = 0; path_index < paths.size(); path_index++) {
            const auto& path = paths[path_index];

            mgp::Path result_path = TranslatePath(graph, graph_view, path);
            mgp::Node current_node = source_node;
            auto costs = ToList(path.costs);

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnNegativeCycle).c_str(), false);
            record.Insert(std::string(shortest_paths::kReturnIndex).c_str(), int64_t(path_index));
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), source_node);
            record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), target_node);
            record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), path.total_cost);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
        }
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}

void DisjointKShortestPaths(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        Source, Sink, K, WeightProp, ScoreProp, DefaultWeight, DefaultScore, CullPerRound, CullAscending
    };

    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        const auto source_node = arguments[ArgIdx::Source].ValueNode();
        const auto sink_node = arguments[ArgIdx::Sink].ValueNode();
        const auto K = arguments[ArgIdx::K].ValueInt();
        const auto weight_property = arguments[ArgIdx::WeightProp].ValueString();
        const auto score_property = std::string(arguments[ArgIdx::ScoreProp].ValueString());
        const auto default_weight = arguments[ArgIdx::DefaultWeight].ValueDouble();
        const auto default_score = arguments[ArgIdx::DefaultScore].ValueDouble();
        const auto cull_per_round = arguments[ArgIdx::CullPerRound].ValueInt();
        const auto cull_ascending = arguments[ArgIdx::CullAscending].ValueBool();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.data() : nullptr;

        if (K < 0) {
            throw mgp::ValueException("K cannot be negative");
        }
        if (!score_property.empty() && cull_per_round < 1) {
            throw mgp::ValueException("cull_per_round must be at least 1");
        }

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );

        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        const auto source_mgid = source_node.Id().AsUint();
        const auto sink_mgid = sink_node.Id().AsUint();

        const auto maybe_source_id = GetInnerNodeIdOpt(graph_view, source_mgid);
        const auto maybe_sink_id = GetInnerNodeIdOpt(graph_view, sink_mgid);
        if (!maybe_source_id || !maybe_sink_id) {
            // Either source or sink is not in subgraph, so there can be no paths between them.
            return;
        }
        const auto source_id = maybe_source_id.value();
        const auto sink_id = maybe_sink_id.value();

        auto abort_func = [&graph] () { graph.CheckMustAbort(); };
        std::vector<shortest_paths::Path<>> paths;
        if (!score_property.empty()) {
            const auto scores = GetScores(graph, graph_view, score_property, default_score);
            paths = shortest_paths::PartialDisjointKShortestPaths(graph_view, source_id, sink_id, K, scores, cull_per_round, cull_ascending, abort_func);
        } else {
            paths = shortest_paths::DisjointKShortestPaths(graph_view, source_id, sink_id, K, abort_func);
        }

        for (std::size_t path_index = 0; path_index < paths.size(); path_index++) {
            const auto& path = paths[path_index];

            mgp::Path result_path = TranslatePath(graph, graph_view, path);
            mgp::Node current_node = source_node;
            auto costs = ToList(path.costs);

            auto record = record_factory.NewRecord();
            record.Insert(std::string(shortest_paths::kReturnIndex).c_str(), int64_t(path_index));
            record.Insert(std::string(shortest_paths::kReturnSourceNode).c_str(), source_node);
            record.Insert(std::string(shortest_paths::kReturnTargetNode).c_str(), sink_node);
            record.Insert(std::string(shortest_paths::kReturnTotalCost).c_str(), path.total_cost);
            record.Insert(std::string(shortest_paths::kReturnCosts).c_str(), costs);
            record.Insert(std::string(shortest_paths::kReturnPath).c_str(), result_path);
        }
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}

} // namespace

extern "C" int mgp_init_module(struct mgp_module *module, struct mgp_memory *memory) {
    try {
        mgp::MemoryDispatcherGuard mem_guard(memory);

        auto return_costs_type = std::make_pair(mgp::Type::List, mgp::Type::Double);
        mgp::Value empty_list(mgp::List{});

        mgp::AddProcedure(
            YensKShortestPaths, shortest_paths::kProcedureYens, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentSourceNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentTargetNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentK, mgp::Type::Int, 1L),
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentThreads, mgp::Type::Int, 0L),
            },
            {
                mgp::Return(shortest_paths::kReturnIndex, mgp::Type::Int),
                mgp::Return(shortest_paths::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTargetNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTotalCost, mgp::Type::Double),
                mgp::Return(shortest_paths::kReturnCosts, return_costs_type),
                mgp::Return(shortest_paths::kReturnPath, mgp::Type::Path)
            },
            module, memory
        );

        mgp::AddProcedure(
            BellmanFordProcedure, shortest_paths::kProcedureBellmanFord, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentSourceNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentTargetNodes, {mgp::Type::List, mgp::Type::Node}, empty_list),
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0)
            },
            {
                mgp::Return(shortest_paths::kReturnNegativeCycle, mgp::Type::Bool),
                mgp::Return(shortest_paths::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTargetNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTotalCost, mgp::Type::Double),
                mgp::Return(shortest_paths::kReturnCosts, return_costs_type),
                mgp::Return(shortest_paths::kReturnPath, mgp::Type::Path)
            },
            module, memory
        );

        mgp::AddProcedure(
            IterativeBellmanFordProcedure, shortest_paths::kProcedureIterativeBellmanFordTargeted, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentSourceNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentTargetNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentRelationshipScoreProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentDefaultScore, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentCullAscending, mgp::Type::Bool, true)
            },
            {
                mgp::Return(shortest_paths::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTargetNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTotalCost, mgp::Type::Double),
                mgp::Return(shortest_paths::kReturnCosts, return_costs_type),
                mgp::Return(shortest_paths::kReturnNumEdgesRemoved, mgp::Type::Int),
                mgp::Return(shortest_paths::kReturnPath, mgp::Type::Path)
            },
            module, memory
        );

        mgp::AddProcedure(
            IBF_SSSP, shortest_paths::kProcedureIterativeBellmanFordSssp, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentSourceNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentTargetNodes, {mgp::Type::List, mgp::Type::Node}, empty_list),
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentRelationshipScoreProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentDefaultScore, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentCullAscending, mgp::Type::Bool, true)
            },
            {
                mgp::Return(shortest_paths::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTargetNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTotalCost, mgp::Type::Double),
                mgp::Return(shortest_paths::kReturnCosts, return_costs_type),
                mgp::Return(shortest_paths::kReturnPath, mgp::Type::Path)
            },
            module, memory
        );

        mgp::AddProcedure(
            IBF_SSSP_Subgraph, shortest_paths::kProcedureIterativeBellmanFordSubgraph, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentSourceNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentRelationshipScoreProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentDefaultScore, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentCullAscending, mgp::Type::Bool, true)
            },
            {
                mgp::Return(shortest_paths::kReturnNodes, {mgp::Type::List, mgp::Type::Node}),
                mgp::Return(shortest_paths::kReturnEdges, {mgp::Type::List, mgp::Type::Relationship}),
                mgp::Return(shortest_paths::kReturnNumEdgesRemoved, mgp::Type::Int),
                mgp::Return(shortest_paths::kReturnRemovedEdges, {mgp::Type::List, mgp::Type::Relationship}),
            },
            module, memory
        );

        mgp::AddProcedure(
            Johnsons_Subgraph, shortest_paths::kProcedureJohnsonsSubgraph, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentRelationshipScoreProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentDefaultScore, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentCullAscending, mgp::Type::Bool, true),
                mgp::Parameter(shortest_paths::kArgumentThreads, mgp::Type::Int, 0L)
            }, {
                mgp::Return(shortest_paths::kReturnNegativeCycle, mgp::Type::Bool),
                mgp::Return(shortest_paths::kReturnNodes, {mgp::Type::List, mgp::Type::Node}),
                mgp::Return(shortest_paths::kReturnEdges, {mgp::Type::List, mgp::Type::Relationship}),
                mgp::Return(shortest_paths::kReturnNodeWeights, {mgp::Type::List, mgp::Type::Double}),
                mgp::Return(shortest_paths::kReturnEdgeWeights, {mgp::Type::List, mgp::Type::Double}),
                mgp::Return(shortest_paths::kReturnNumEdgesRemoved, mgp::Type::Int),
                mgp::Return(shortest_paths::kReturnRemovedEdges, {mgp::Type::List, mgp::Type::Relationship}),
            },
            module, memory
        );

        mgp::AddProcedure(
            Johnsons_Source_Subgraphs, shortest_paths::kProcedureJohnsonsSourceSubgraphs, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentSourceNodes, {mgp::Type::List, mgp::Type::Node}, empty_list),
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentRelationshipScoreProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentDefaultScore, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentCullAscending, mgp::Type::Bool, true),
                mgp::Parameter(shortest_paths::kArgumentThreads, mgp::Type::Int, 0L)
            }, {
                mgp::Return(shortest_paths::kReturnNegativeCycle, mgp::Type::Bool),
                mgp::Return(shortest_paths::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnNodes, {mgp::Type::List, mgp::Type::Node}),
                mgp::Return(shortest_paths::kReturnEdges, {mgp::Type::List, mgp::Type::Relationship}),
                mgp::Return(shortest_paths::kReturnCosts, {mgp::Type::List, mgp::Type::Double})
            },
            module, memory
        );

        mgp::AddProcedure(
            Johnsons_Paths, shortest_paths::kProcedureJohnsonsPaths, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentSourceNodes, {mgp::Type::List, mgp::Type::Node}, empty_list),
                mgp::Parameter(shortest_paths::kArgumentTargetNodes, {mgp::Type::List, mgp::Type::Node}, empty_list),
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentRelationshipScoreProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentDefaultScore, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentCullAscending, mgp::Type::Bool, true),
                mgp::Parameter(shortest_paths::kArgumentThreads, mgp::Type::Int, 0L)
            }, {
                mgp::Return(shortest_paths::kReturnNegativeCycle, mgp::Type::Bool),
                mgp::Return(shortest_paths::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTargetNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTotalCost, mgp::Type::Double),
                mgp::Return(shortest_paths::kReturnCosts, return_costs_type),
                mgp::Return(shortest_paths::kReturnPath, mgp::Type::Path)
            },
            module, memory
        );

        mgp::AddProcedure(
            Johnsons_K_Shortest, shortest_paths::kProcedureJohnsonsKShortest, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentSourceNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentTargetNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentK, mgp::Type::Int, 1L),
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentRelationshipScoreProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentDefaultScore, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentCullAscending, mgp::Type::Bool, true),
                mgp::Parameter(shortest_paths::kArgumentThreads, mgp::Type::Int, 0L)
            }, {
                mgp::Return(shortest_paths::kReturnNegativeCycle, mgp::Type::Bool),
                mgp::Return(shortest_paths::kReturnIndex, mgp::Type::Int),
                mgp::Return(shortest_paths::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTargetNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTotalCost, mgp::Type::Double),
                mgp::Return(shortest_paths::kReturnCosts, return_costs_type),
                mgp::Return(shortest_paths::kReturnPath, mgp::Type::Path)
            },
            module, memory
        );

        mgp::AddProcedure(
            Johnsons_Disjoint_K_Shortest, shortest_paths::kProcedureJohnsonsDisjointKShortest, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentSourceNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentTargetNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentK, mgp::Type::Int, 0L),
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentRelationshipScoreProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentDefaultScore, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentCullAscending, mgp::Type::Bool, true)
            }, {
                mgp::Return(shortest_paths::kReturnNegativeCycle, mgp::Type::Bool),
                mgp::Return(shortest_paths::kReturnIndex, mgp::Type::Int),
                mgp::Return(shortest_paths::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTargetNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTotalCost, mgp::Type::Double),
                mgp::Return(shortest_paths::kReturnCosts, return_costs_type),
                mgp::Return(shortest_paths::kReturnPath, mgp::Type::Path)
            },
            module, memory
        );

        mgp::AddProcedure(
            DisjointKShortestPaths, shortest_paths::kProcedureDisjointKShortest, mgp::ProcedureType::Read,
            {
                mgp::Parameter(shortest_paths::kArgumentSourceNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentTargetNode, mgp::Type::Node),
                mgp::Parameter(shortest_paths::kArgumentK, mgp::Type::Int, 0L),
                mgp::Parameter(shortest_paths::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentRelationshipScoreProperty, mgp::Type::String, ""),
                mgp::Parameter(shortest_paths::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentDefaultScore, mgp::Type::Double, 1.0),
                mgp::Parameter(shortest_paths::kArgumentCullPerRound, mgp::Type::Int, 1L),
                mgp::Parameter(shortest_paths::kArgumentCullAscending, mgp::Type::Bool, true)
            },
            {
                mgp::Return(shortest_paths::kReturnIndex, mgp::Type::Int),
                mgp::Return(shortest_paths::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTargetNode, mgp::Type::Node),
                mgp::Return(shortest_paths::kReturnTotalCost, mgp::Type::Double),
                mgp::Return(shortest_paths::kReturnCosts, return_costs_type),
                mgp::Return(shortest_paths::kReturnPath, mgp::Type::Path)
            },
            module, memory
        );

    } catch (const std::exception& e) {
        return 1;
    }

    return 0;
}

extern "C" int mgp_shutdown_module() { return 0; }