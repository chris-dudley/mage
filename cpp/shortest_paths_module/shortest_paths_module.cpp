#include <mgp.hpp>
#include <mg_utils.hpp>

#include <stdexcept>
#include <string>

#include "algorithm/shortest_path.hpp"
#include "algorithm/yens.hpp"
#include "algorithm/procedures.hpp"
#include "algorithm/bellman_ford.hpp"
#include "algorithm/iterative_bf.hpp"

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

void YensKShortestPaths(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    // Indicies in the arguments list for parameters
    enum ArgIdx : size_t {
        Source, Sink, K, WeightProp, DefaultWeight
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

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.data() : nullptr;

        if (K < 0) {
            throw mgp::ValueException("K cannot be negative");
        }
        if (K == 0) {
            // Why did you even call this procedure?
            return;
        }

        shortest_paths::ShortestPathFunc sp_func = shortest_paths::Dijkstra;

        auto graph_view_ptr = mg_utility::GetGraphView(
            memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
            weighted, weight_prop_cstr, default_weight
        );

        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        const auto source_mgid = source_node.Id().AsUint();
        const auto sink_mgid = sink_node.Id().AsUint();

        const auto source_id = GetInnerNodeId(graph_view, source_mgid);
        const auto sink_id = GetInnerNodeId(graph_view, sink_mgid);

        auto abort_func = [&graph] () { graph.CheckMustAbort(); };
        const auto paths = shortest_paths::KShortestPaths(graph_view, source_id, sink_id, K, sp_func, abort_func);

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
        const auto source_id = GetInnerNodeId(graph_view, source_mgid);

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
        const auto source_id = GetInnerNodeId(graph_view, source_mgid);
        const auto target_mgid = target_node.Id().AsUint();
        const auto target_id = GetInnerNodeId(graph_view, target_mgid);
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
        const auto source_id = GetInnerNodeId(graph_view, source_mgid);
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
        const auto source_id = GetInnerNodeId(graph_view, source_mgid);
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
    } catch (const std::exception& e) {
        return 1;
    }

    return 0;
}

extern "C" int mgp_shutdown_module() { return 0; }