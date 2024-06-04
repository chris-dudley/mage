#pragma once

#include <vector>

#include <mg_graph.hpp>
#include "Highs.h"

#include "path.hpp"
#include "polyedge.hpp"

namespace shortest_paths {

struct OptimizationResult {
    std::vector<double> input_ratios;
    std::vector<double> input_amounts;
    double expected_output;
    bool success;
};

// TODO: Placeholder interface
template <typename TSize = std::uint64_t>
OptimizationResult OptimizeFlows(const mg_graph::GraphView<TSize>& graph, const std::vector<Path<TSize>>& paths) {
    OptimizationResult result {
        std::vector<double>(paths.size(), 0.0),
        std::vector<double>(paths.size(), 0.0),
        0.0,
        false
    };

    // TODO: Actually implement
    result.success = true;

    return result;
}

} // namespace shortest_paths