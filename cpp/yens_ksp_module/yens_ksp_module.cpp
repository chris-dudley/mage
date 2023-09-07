#include <mgp.hpp>
#include <mg_utils.hpp>

#include <stdexcept>
#include <string>

#include "algorithm/shortest_path.hpp"
#include "algorithm/yens.hpp"

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

// Add wrappers for GetMemgraph(Node|Edge)Id just in case.
uint64_t GetMemgraphNodeId(const mg_graph::GraphView<>& graph, uint64_t inner_id) {
    return graph.GetMemgraphNodeId(inner_id);
}
uint64_t GetMemgraphEdgeId(const mg_graph::GraphView<>& graph, uint64_t inner_id) {
    return graph.GetMemgraphEdgeId(inner_id);
}

void YensKShortestPaths(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory, bool subgraph) {
    mgp::MemoryDispatcherGuard mem_guard(memory);
    mgp::Graph graph{memgraph_graph};

    const auto arguments = mgp::List(args);
    const auto record_factory = mgp::RecordFactory(result);
    try {
        size_t args_idx = 0;
        mgp_list *subgraph_nodes = nullptr, *subgraph_edges = nullptr;
        if (subgraph) {
            // Can't use wrappers here, as we need the underlying mgp_list* values to pass to
            // GetSubgraphView
            subgraph_nodes = mgp::value_get_list(mgp::list_at(args, args_idx++));
            subgraph_edges = mgp::value_get_list(mgp::list_at(args, args_idx++));
        }

        const auto source_node = arguments[args_idx++].ValueNode();
        const auto sink_node = arguments[args_idx++].ValueNode();
        const auto K = arguments[args_idx++].ValueInt();
        const auto weight_property = arguments[args_idx++].ValueString();
        const auto default_weight = arguments[args_idx++].ValueDouble();

        bool weighted = ! weight_property.empty();
        const char * weight_prop_cstr = weighted ? weight_property.data() : nullptr;

        if (K < 0) {
            throw mgp::ValueException("K cannot be negative");
        }
        if (K == 0) {
            // Why did you even call this procedure?
            return;
        }

        yens_alg::ShortestPathFunc sp_func = yens_alg::Dijkstra;

        auto graph_view_ptr = subgraph ?
            mg_utility::GetSubgraphView(
                memgraph_graph, result, memory, subgraph_nodes, subgraph_edges, mg_graph::GraphType::kDirectedGraph,
                weighted, weight_prop_cstr, default_weight)
            :
            mg_utility::GetGraphView(
                memgraph_graph, result, memory, mg_graph::GraphType::kDirectedGraph,
                weighted, weight_prop_cstr, default_weight);

        const mg_graph::GraphView<>& graph_view = *graph_view_ptr;

        const auto source_mgid = source_node.Id().AsUint();
        const auto sink_mgid = sink_node.Id().AsUint();

        const auto source_id = GetInnerNodeId(graph_view, source_mgid);
        const auto sink_id = GetInnerNodeId(graph_view, sink_mgid);

        auto abort_func = [&graph] () { graph.CheckMustAbort(); };
        const auto paths = yens_alg::KShortestPaths(graph_view, source_id, sink_id, K, sp_func, abort_func);

        for (std::size_t path_index = 0; path_index < paths.size(); path_index++) {
            const auto& path = paths[path_index];

            mgp::Path result_path(source_node);
            mgp::Node current_node = source_node;
            mgp::List costs(path.weights.size()); // accumulated cost at each node

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
            for (auto cost : path.weights) {
                costs.Append(mgp::Value(cost));
            }

            auto record = record_factory.NewRecord();
            record.Insert(std::string(yens_alg::kReturnIndex).c_str(), int64_t(path_index));
            record.Insert(std::string(yens_alg::kReturnSourceNode).c_str(), source_node);
            record.Insert(std::string(yens_alg::kReturnTargetNode).c_str(), sink_node);
            record.Insert(std::string(yens_alg::kReturnTotalCost).c_str(), path.total_weight);
            record.Insert(std::string(yens_alg::kReturnCosts).c_str(), costs);
            record.Insert(std::string(yens_alg::kReturnPath).c_str(), result_path);
        }
    } catch (const std::exception& e) {
        record_factory.SetErrorMessage(e.what());
        return;
    }
}


void YensGraph(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    YensKShortestPaths(args, memgraph_graph, result, memory, false);
}

void YensSubgraph(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
    YensKShortestPaths(args, memgraph_graph, result, memory, true);
}


} // namespace

extern "C" int mgp_init_module(struct mgp_module *module, struct mgp_memory *memory) {
    try {
        mgp::MemoryDispatcherGuard mem_guard(memory);

        auto return_costs_type = std::make_pair(mgp::Type::List, mgp::Type::Double);
        auto subgraph_nodes_type = std::make_pair(mgp::Type::List, mgp::Type::Node);
        auto subgraph_edges_type = std::make_pair(mgp::Type::List, mgp::Type::Relationship);

        mgp::AddProcedure(
            YensGraph, yens_alg::kProcedureKShortestPaths, mgp::ProcedureType::Read,
            {
                mgp::Parameter(yens_alg::kArgumentSourceNode, mgp::Type::Node),
                mgp::Parameter(yens_alg::kArgumentTargetNode, mgp::Type::Node),
                mgp::Parameter(yens_alg::kArgumentK, mgp::Type::Int, 1L),
                mgp::Parameter(yens_alg::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(yens_alg::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
            },
            {
                mgp::Return(yens_alg::kReturnIndex, mgp::Type::Int),
                mgp::Return(yens_alg::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(yens_alg::kReturnTargetNode, mgp::Type::Node),
                mgp::Return(yens_alg::kReturnTotalCost, mgp::Type::Double),
                mgp::Return(yens_alg::kReturnCosts, return_costs_type),
                mgp::Return(yens_alg::kReturnPath, mgp::Type::Path)
            },
            module, memory
        );

        mgp::AddProcedure(
            YensSubgraph, yens_alg::kProcedureSubgraphKShortestPaths, mgp::ProcedureType::Read,
            {
                mgp::Parameter(yens_alg::kArgumentSubgraphNodes, subgraph_nodes_type),
                mgp::Parameter(yens_alg::kArgumentSubgraphEdges, subgraph_edges_type),
                mgp::Parameter(yens_alg::kArgumentSourceNode, mgp::Type::Node),
                mgp::Parameter(yens_alg::kArgumentTargetNode, mgp::Type::Node),
                mgp::Parameter(yens_alg::kArgumentK, mgp::Type::Int, 1L),
                mgp::Parameter(yens_alg::kArgumentRelationshipWeightProperty, mgp::Type::String, ""),
                mgp::Parameter(yens_alg::kArgumentDefaultWeight, mgp::Type::Double, 1.0),
            },
            {
                mgp::Return(yens_alg::kReturnIndex, mgp::Type::Int),
                mgp::Return(yens_alg::kReturnSourceNode, mgp::Type::Node),
                mgp::Return(yens_alg::kReturnTargetNode, mgp::Type::Node),
                mgp::Return(yens_alg::kReturnTotalCost, mgp::Type::Double),
                mgp::Return(yens_alg::kReturnCosts, return_costs_type),
                mgp::Return(yens_alg::kReturnPath, mgp::Type::Path)
            },
            module, memory
        );

    } catch (const std::exception& e) {
        return 1;
    }

    return 0;
}

extern "C" int mgp_shutdown_module() { return 0; }