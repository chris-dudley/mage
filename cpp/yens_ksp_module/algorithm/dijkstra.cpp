#include "shortest_path.hpp"

#include <ranges>
#include <limits>
#include <queue>
#include <functional>
#include <unordered_map>
#include <stdexcept>

namespace yens_alg {

using std::uint64_t;

constexpr const uint64_t INVALID_NODE = std::numeric_limits<uint64_t>::max();
constexpr const uint64_t INVALID_EDGE = std::numeric_limits<uint64_t>::max();
constexpr const double INFINITY = std::numeric_limits<double>::infinity();


using NodeAndDistance = std::pair<uint64_t, double>;

bool NodeAndDistanceGreater(const NodeAndDistance &a, const NodeAndDistance &b) {
    return a.second > b.second;
}

Path<> Dijkstra(
    const mg_graph::GraphView<> &graph, uint64_t source_id, uint64_t sink_id,
    const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes
) {
    Path<> result{source_id};
    if (source_id == sink_id) {
        return result;
    }

    const auto &verticies = graph.Nodes();

    // `shortest_distance[i]` holds the shortest distance from `source` to `i`.
    std::vector<double> shortest_distance(verticies.size(), INFINITY);
    // `added[i]` will be true if vertex `i` is included in the shortest path tree or
    // shortest distance from `source` to `i` is finalized.
    std::vector<bool> added(verticies.size(), false);
    // `parent[i]` is the ID of the parent node of vertex `i` in the shortest path tree.
    std::vector<uint64_t> parent(verticies.size(), INVALID_NODE);
    // `edge_in[i]` is the ID of the edge used to reach node `i` in the shortest path tree.
    std::vector<uint64_t> edge_in(verticies.size(), INVALID_EDGE);

    // Distance of source vertex to itself is always 0
    shortest_distance[source_id] = 0.0;

    // node_queue contains a queue of nodes to visit, sorted by the total cost to
    // reach that node. The top of the queue is the node with the minimum distance.
    std::priority_queue<
        // std::greater makes it a min queue
        NodeAndDistance, std::vector<NodeAndDistance>, std::greater<NodeAndDistance>
    > node_queue;

    // Seed queue with source node
    node_queue.emplace(source_id, 0.0);

    while (!node_queue.empty()) {
        auto current_node = node_queue.top().first;
        node_queue.pop();

        if (current_node == sink_id) {
            // Found the target, we can stop now
            break;
        }

        if (added[current_node]) {
            // Don't revisit nodes.
            continue;
        }
        added[current_node] = true;

        for (auto neighbor : graph.OutNeighbours(current_node)) {
            if (ignored_edges.contains(neighbor.edge_id) || ignored_nodes.contains(neighbor.node_id)) {
                continue;
            }
            double edge_weight = graph.IsWeighted() ? graph.GetWeight(neighbor.edge_id) : 1.0;
            double neighbor_weight = shortest_distance[current_node] + edge_weight;

            if (neighbor_weight < shortest_distance[neighbor.node_id]) {
                shortest_distance[neighbor.node_id] = neighbor_weight;

                parent[neighbor.node_id] = current_node;
                edge_in[neighbor.node_id] = neighbor.edge_id;
                node_queue.emplace(neighbor.node_id, neighbor_weight);
            }
        }
    }

    if(parent[sink_id] == INVALID_EDGE) {
        // Return empty path
        return result;
    }

    std::vector<uint64_t> edges_stack;
    for (uint64_t current_node = sink_id; current_node != source_id; current_node = parent[current_node]) {
        edges_stack.push_back(edge_in[current_node]);
    }

    for (auto edge_id : std::views::reverse(edges_stack)) {
        result.add_edge(graph.GetEdge(edge_id), graph.IsWeighted() ? graph.GetWeight(edge_id) : 1.0);
    }
    return result;
}

} // namespace yens_alg