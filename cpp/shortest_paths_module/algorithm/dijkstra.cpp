#include "shortest_path.hpp"
#include "dijkstra.hpp"

#include <ranges>
#include <limits>
#include <queue>
#include <functional>
#include <unordered_map>
#include <stdexcept>

namespace shortest_paths {

void CheckAbortNoop() { }

Path<> Dijkstra(
    const mg_graph::GraphView<> &graph, uint64_t source_id, uint64_t sink_id,
    const EdgeIdSet& ignored_edges, const NodeIdSet& ignored_nodes,
    const CheckAbortFunc& check_abort
) {
    DijkstraPathfinder<> pathfinder;
    pathfinder.search(graph, source_id, sink_id, ignored_edges, ignored_nodes, check_abort);
    return pathfinder.path_to(sink_id);
}

} // namespace shortest_paths