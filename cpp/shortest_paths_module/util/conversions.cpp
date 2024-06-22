#include "conversions.hpp"

#include <cmath>
#include <algorithm>

#include <fmt/core.h>

namespace shortest_paths {
namespace util {

EdgeNetwork<uint64_t> EdgeNetworkFromMemgraph(
    const mgp::Graph& graph, std::vector<mgp::Path> paths,
    const std::string& x_coord_prop_name, const std::string& y_coord_prop_name
) {
    std::vector<std::vector<PolyEdge<uint64_t>>> polyedge_paths;

    for (size_t path_idx = 0; path_idx < paths.size(); path_idx++) {
        const auto& path = paths[path_idx];
        std::vector<PolyEdge<uint64_t>> cur_path;
        for (size_t edge_idx = 0; edge_idx < path.Length(); edge_idx++) {
            auto edge = path.GetRelationshipAt(edge_idx);

            auto edge_id = edge.Id().AsUint();

            std::vector<double> x_coords;
            std::vector<double> y_coords;

            // Validate and extract x coordinates
            auto x_coords_value = edge.GetProperty(x_coord_prop_name);
            if (!x_coords_value.IsList()) {
                throw std::invalid_argument(fmt::format(
                    "path {}, edge {} (id: {}) property '{}' is not a list",
                    path_idx, edge_idx, edge_id, x_coord_prop_name
                ));
            }
            auto x_coords_list = x_coords_value.ValueList();
            for (auto x_value : x_coords_list) {
                if (x_value.IsNumeric()) {
                    x_coords.push_back(x_value.ValueNumeric());
                } else if(x_value.IsNull()) {
                    x_coords.push_back(NAN);
                } else {
                    throw std::invalid_argument(fmt::format(
                        "path {}, edge {} (id: {}) property '{}' contains non-numeric value at index {}",
                        path_idx, edge_idx, edge_id, x_coord_prop_name, x_coords.size()
                    ));
                }
            }

            // Validate and extract y coordinates
            auto y_coords_value = edge.GetProperty(y_coord_prop_name);
            if (!y_coords_value.IsList()) {
                throw std::invalid_argument(fmt::format(
                    "path {}, edge {} (id: {}) property '{}' is not a list",
                    path_idx, edge_idx, edge_id, y_coord_prop_name
                ));
            }
            for (auto y_value : x_coords_list) {
                if (y_value.IsNumeric()) {
                    y_coords.push_back(y_value.ValueNumeric());
                } else if(y_value.IsNull()) {
                    y_coords.push_back(NAN);
                } else {
                    throw std::invalid_argument(fmt::format(
                        "path {}, edge {} (id: {}) property '{}' contains non-numeric value at index {}",
                        path_idx, edge_idx, edge_id, y_coord_prop_name, y_coords.size()
                    ));
                }
            }

            // Find out how long a prefix of finite x and y coordinate pairs we have. We'll automatically strip of
            // suffixes with non-finite values as long as we have a prefix of at least 2 finite pairs.
            size_t finite_prefix_len;
            for (size_t i = 0; i < std::min(x_coords.size(), y_coords.size()); i++) {
                if (std::isfinite(x_coords[i]) && std::isfinite(y_coords[i])) {
                    finite_prefix_len++;
                } else {
                    break;
                }
            }

            if (finite_prefix_len < 2) {
                throw std::invalid_argument(fmt::format(
                    "path {}, edge {} (id: {}) has only prefix of {} finite x and y coordinate pairs, need at least 2",
                        path_idx, edge_idx, edge_id, finite_prefix_len
                ));
            }
            if (x_coords.size() > finite_prefix_len) {
                x_coords.resize(finite_prefix_len);
            }
            if (y_coords.size() > finite_prefix_len) {
                y_coords.resize(finite_prefix_len);
            }

            cur_path.emplace_back(edge_id, x_coords, y_coords);
        }

        polyedge_paths.emplace_back(std::move(cur_path));
    }

    return EdgeNetwork<uint64_t>(polyedge_paths);
}

} // namespace util
} // namespace shortest_paths