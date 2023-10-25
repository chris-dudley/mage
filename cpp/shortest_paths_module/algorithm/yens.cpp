#include "yens.hpp"
#include "dijkstra.hpp"

#include <algorithm>
#include <queue>
#include <functional>

namespace shortest_paths {

std::vector<Path<>> KShortestPaths(
    const mg_graph::GraphView<> &graph, std::uint64_t source_id, std::uint64_t sink_id,
    std::uint64_t K, const CheckAbortFunc& check_abort, int threads
) {
    YensPathfinder<> pathfinder;
    return pathfinder.search(graph, source_id, sink_id, K, check_abort, threads);
}

} // namespace shortest_paths