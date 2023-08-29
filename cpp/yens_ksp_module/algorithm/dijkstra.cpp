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

struct NodeAndDistance {
    uint64_t NodeId;
    double Weight;
};

constexpr double operator<=>( const NodeAndDistance& lhs, const NodeAndDistance& rhs ) {
    return lhs.Weight - rhs.Weight;
}

bool NodeAndDistanceGreater(const NodeAndDistance &a, const NodeAndDistance &b) {
    return a.Weight > b.Weight;
}

template<typename K, typename V>
const V& get_default(const std::unordered_map<K, V>& map, const K& key, const V& default_value) {
    auto iter = map.find(key);
    if (iter != map.end()) {
        return iter->second;
    }
    return default_value;
}

Path<> Dijkstra(
    const mg_graph::GraphView<> &graph, uint64_t source_id, uint64_t sink_id,
    const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes
) {
    Path<> result{source_id};
    if (source_id == sink_id) {
        return result;
    }

    // Distance from start node to this one
    std::unordered_map<uint64_t, double> distances;
    // ID of node used to get to this one
    std::unordered_map<uint64_t, uint64_t> previous;
    // ID of edge used to get to this one
    std::unordered_map<uint64_t, uint64_t> edges;
    std::priority_queue<
        // std::greater makes it a min queue
        NodeAndDistance, std::vector<NodeAndDistance>, std::greater<NodeAndDistance>
    > node_queue;

    distances[source_id] = 0.0;

    // Seed queue with source node
    node_queue.emplace(source_id, 0.0);

    while (!node_queue.empty()) {
        auto current = node_queue.top();
        node_queue.pop();

        if (current.NodeId == sink_id) {
            // Found the target, we can stop now
            break;
        }

        double current_distance = get_default(distances, current.NodeId, INFINITY);

        if (current.Weight > current_distance) {
            // Node+Weight was added before another, better path was found, ignore.
            continue;
        }

        for (auto neighbor : graph.OutNeighbours(current.NodeId)) {
            if (ignored_edges.contains(neighbor.edge_id) || ignored_nodes.contains(neighbor.node_id)) {
                continue;
            }
            double edge_weight = graph.GetWeight(neighbor.edge_id);
            double neighbor_weight = current_distance + edge_weight;

            if (neighbor_weight < get_default(distances, neighbor.node_id, INFINITY)) {
                distances[neighbor.node_id] = neighbor_weight;
                previous[neighbor.node_id] = current.NodeId;
                edges[neighbor.node_id] = neighbor.edge_id;
                node_queue.emplace(neighbor.node_id, neighbor_weight);
            }
        }
    }

    if(!previous.contains(sink_id)) {
        // Return empty path
        return result;
    }

    uint64_t current_node = sink_id;
    std::vector<uint64_t> edges_stack;
    while (current_node != source_id) {
        edges_stack.push_back(edges.at(current_node));
        current_node = previous.at(current_node);
    }

    for (auto edge_id : std::views::reverse(edges_stack)) {
        result.add_edge(graph.GetEdge(edge_id), graph.GetWeight(edge_id));
    }
    return result;
}

} // namespace yens_alg