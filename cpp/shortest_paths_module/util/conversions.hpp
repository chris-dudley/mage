#pragma once

#include <mgp.hpp>

#include "../algorithm/polyedge.hpp"

namespace shortest_paths {
namespace util {

EdgeNetwork<uint64_t> EdgeNetworkFromMemgraph(const mgp::Graph &graph, std::vector<mgp::Path> paths,
                                              const std::string &x_coord_prop_name,
                                              const std::string &y_coord_prop_name, bool fit_2point_quadratic = false);

}  // namespace util
}  // namespace shortest_paths
