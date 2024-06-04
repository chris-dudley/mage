#include <map>
#include <set>
#include <unordered_set>
#include <list>
#include <forward_list>

#include <gtest/gtest.h>
#include <mg_generate.hpp>
#include <mg_graph.hpp>
#include <mg_test_utils.hpp>

#include "algorithm/shortest_path.hpp"
#include "algorithm/dijkstra.hpp"
#include "algorithm/yens.hpp"
#include "algorithm/bellman_ford.hpp"
#include "algorithm/iterative_bf.hpp"
#include "algorithm/johnsons.hpp"
#include "algorithm/disjoint.hpp"
#include "algorithm/successive_shortest_paths.hpp"
#include "algorithm/polyedge.hpp"
#include "algorithm/optimize.hpp"

void check_abort_noop() {}

bool CheckPaths(const std::vector<shortest_paths::Path<>>& paths, const std::vector<shortest_paths::Path<>>& expected) {
    std::multimap<double, size_t> cost_to_path_index;
    for (size_t i = 0; i < expected.size(); i++) {
        cost_to_path_index.emplace(expected[i].total_cost, i);
    }

    // Paths of equal length could be returned in any order and still be valid
    for (const auto& path : paths) {
        auto range = cost_to_path_index.equal_range(path.total_cost);
        if (range.first == range.second) {
            // No match
            return false;
        }

        // Search for a matching path in the expected list
        bool found_match = false;
        for (auto iter = range.first; iter != range.second; iter++) {
            const auto& other = expected[iter->second];

            if (path == other) {
                // Found match, remove it so that we can't match it again.
                found_match = true;
                cost_to_path_index.erase(iter);
                break;
            }
        }

        if (!found_match) {
            return false;
        }
    }

    return true;
}

TEST(ShortestPaths, DijkstraEmptyGraph) {
    auto G = mg_generate::BuildGraph(0, {});

    auto path = shortest_paths::Dijkstra(*G, 0, 0, {}, {}, check_abort_noop);

    ASSERT_TRUE(path.empty());
}

TEST(ShortestPaths, DijkstraSingleNode) {
    auto G = mg_generate::BuildGraph(1, {});

    auto path = shortest_paths::Dijkstra(*G, 0, 0, {}, {}, check_abort_noop);

    ASSERT_TRUE(path.empty());
}

TEST(ShortestPaths, DijkstraDisconnectedNodes) {
    auto G = mg_generate::BuildGraph(10, {});

    auto path = shortest_paths::Dijkstra(*G, 0, 9, {}, {}, check_abort_noop);

    ASSERT_TRUE(path.empty());
}

/*
 *                ┌───┐                          ┌───┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, DijkstraSmallAcyclicGraph) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /*0*/ {{0, 1}, 50.0}, /*1*/ {{0, 2}, 50.0}, /*2*/ {{0, 3}, 100.0},
            /*3*/ {{1, 3}, 40.0},
            /*4*/ {{2, 3}, 40.0}, /*5*/ {{2, 4}, 80.0},
            /*6*/ {{3, 4}, 30.0}, /*7*/ {{3, 5}, 80.0},
            /*8*/ {{4, 5}, 40.0}
        },
        mg_graph::GraphType::kDirectedGraph
    );
    shortest_paths::Path<> expected_path{
        {0, 1, 3, 4, 5}, // nodes
        {0, 3, 6, 8},    // edges
        {0.0, 50.0, 90.0, 120.0, 160.0}, // costs
        160.0 // total weight
    };

    auto path = shortest_paths::Dijkstra(*G, 0, 5, {}, {}, check_abort_noop);
    ASSERT_EQ(path, expected_path);
}


/*
 *                ┌───┐                          ┌───┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └─▲─┘         │                └▲─┬┘
 *   │              │          40                 │ │
 *   ├──70──────────┘           │                 │ 40
 *   │                          │                 │ │
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     ├──80─────►│ 5 │
 *   │                          │     │          └─▲─┘
 *   │            ┌───┐         │     │            │
 *   └──50───────►│ 1 ├─────────┘     └──90────────┘
 *                └───┘
 */
TEST(ShortestPaths, DijkstraSmallAcyclicGraphParallelEdges) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, 30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{0, 2}, 70.0},
            /*10*/ {{3, 5}, 90.0}
        },
        mg_graph::GraphType::kDirectedGraph
    );
    // Path should still be the same.
    shortest_paths::Path<> expected_path{
        {0, 1, 3, 4, 5}, // nodes
        {0, 3, 6, 8},    // edges
        {0.0, 50.0, 90.0, 120.0, 160.0}, // costs
        160.0 // total weight
    };

    auto path = shortest_paths::Dijkstra(*G, 0, 5, {}, {}, check_abort_noop);
    ASSERT_EQ(path, expected_path);
}

/*
 *                  ┌─────────────10───────────────┐
 *                ┌─▼─┐                          ┌─┴─┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, DijkstraSmallGraphWithCycle) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, 30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{4, 2}, 10.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );
    // Path should still be the same.
    shortest_paths::Path<> expected_path{
        {0, 1, 3, 4, 5}, // nodes
        {0, 3, 6, 8},    // edges
        {0.0, 50.0, 90.0, 120.0, 160.0}, // costs
        160.0 // total weight
    };

    auto path = shortest_paths::Dijkstra(*G, 0, 5, {}, {}, check_abort_noop);
    ASSERT_EQ(path, expected_path);
}

/*
 *                  ┌─────────── -100 ─────────────┐
 *                ┌─▼─┐                          ┌─┴─┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬─-30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, DijkstraSmallGraphWithNegativeCycle) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, -30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{4, 2}, -100.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );
    // The path isn't correct due to the negative edges, but it should at least return
    // something.
    shortest_paths::Path<> expected_path{
        {0, 1, 3, 4, 5}, // nodes
        {0, 3, 6, 8},    // edges
        {0.0, 50.0, 90.0, 60.0, 100.0}, // costs
        100.0 // total weight
    };

    auto path = shortest_paths::Dijkstra(*G, 0, 5, {}, {}, check_abort_noop);
    ASSERT_EQ(path, expected_path);
}

/*
 *             ┌───┐         ┌───┐
 *   ┌─────────┤ 1 ├─────────┤ 3 ├─────────┐
 *   │         └─┬─┘\       /└─┬─┘         │
 *   │           │   \     /   │           │
 *   │           │    \   /    │           │
 * ┌─┴─┐         │     \ /     │         ┌─┴─┐
 * │ 0 │         │      x      │         │ 5 │
 * └─┬─┘         │     / \     │         └─┬─┘
 *   │           │    /   \    │           │
 *   │         ┌─┴─┐ /     \ ┌─┴─┐         │
 *   └─────────┤ 2 ├─────────┤ 4 ├─────────┘
 *             └───┘         └───┘
 */
TEST(ShortestPaths, DijkstraAllNegativeCycles) {
    // This isn't technically a well-formed graph for Dijstra's, but we just want to
    // make sure it finishes.
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /* 0*/ {{0, 1}, -1.0}, /* 1*/ {{0, 2}, -1.0},
            /* 2*/ {{1, 0}, -1.0}, /* 3*/ {{1, 2}, -1.0}, /* 4*/ {{1, 3}, -1.0}, /* 5*/ {{1, 4}, -1.0},
            /* 6*/ {{2, 0}, -1.0}, /* 7*/ {{2, 1}, -1.0}, /* 8*/ {{2, 3}, -1.0}, /* 9*/ {{2, 4}, -1.0},
            /*10*/ {{3, 1}, -1.0}, /*11*/ {{3, 2}, -1.0}, /*12*/ {{3, 4}, -1.0}, /*13*/ {{3, 5}, -1.0},
            /*14*/ {{4, 1}, -1.0}, /*15*/ {{4, 2}, -1.0}, /*16*/ {{4, 3}, -1.0}, /*17*/ {{4, 5}, -1.0},
            /*18*/ {{5, 3}, -1.0}, /*19*/ {{5, 4}, -1.0}
        },
        mg_graph::GraphType::kDirectedGraph
    );
    auto path = shortest_paths::Dijkstra(*G, 0, 5, {}, {}, check_abort_noop);
    ASSERT_FALSE(path.empty());
    ASSERT_EQ(path.size(), 5ULL);
}

TEST(ShortestPaths, DijkstraHugeCyclicGraph) {
    // Create a graph where each vertex is connected to the next 10 verticies, creating a huge
    // interconnected cycle.
    const uint64_t NUM_VERTICIES = 10'000;
    const uint64_t EDGES_PER_VERTEX = 10;
    const uint64_t SHORTEST_PATH_EDGES = (NUM_VERTICIES / EDGES_PER_VERTEX);
    auto G = mg_generate::BuildGraph(NUM_VERTICIES, {}, mg_graph::GraphType::kDirectedGraph);
    for (uint64_t i = 0; i < NUM_VERTICIES; i++) {
        for (uint64_t mod = 1; mod <= EDGES_PER_VERTEX; mod++) {
            uint64_t next_node = (i + mod) % NUM_VERTICIES;
            G->CreateEdge(i, next_node, mg_graph::GraphType::kDirectedGraph);
        }
    }

    // Find a path from 0 to NUM_VERTICIES-1
    auto path = shortest_paths::Dijkstra(*G, 0, NUM_VERTICIES-1, {}, {}, check_abort_noop);
    ASSERT_EQ(path.empty(), false);
    ASSERT_EQ(path.size(), SHORTEST_PATH_EDGES);
    ASSERT_DOUBLE_EQ(path.total_cost, (SHORTEST_PATH_EDGES) * 1.0);
}

TEST(ShortestPaths, YensEmptyGraph) {
    auto G = mg_generate::BuildGraph(0, {});

    auto paths = shortest_paths::KShortestPaths(*G, 0, 0, 1, check_abort_noop);

    ASSERT_TRUE(paths.empty());
}

TEST(ShortestPaths, YensSingleNode) {
    auto G = mg_generate::BuildGraph(1, {});

    auto paths = shortest_paths::KShortestPaths(*G, 0, 0, 1, check_abort_noop);

    ASSERT_TRUE(paths.empty());
}

TEST(ShortestPaths, YensDisconnectedNodes) {
    auto G = mg_generate::BuildGraph(10, {});

    auto paths = shortest_paths::KShortestPaths(*G, 0, 9, 1, check_abort_noop);

    ASSERT_TRUE(paths.empty());
}

/*
 *                ┌───┐                          ┌───┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, YensSmallAcyclicGraph) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /*0*/ {{0, 1}, 50.0}, /*1*/ {{0, 2}, 50.0}, /*2*/ {{0, 3}, 100.0},
            /*3*/ {{1, 3}, 40.0},
            /*4*/ {{2, 3}, 40.0}, /*5*/ {{2, 4}, 80.0},
            /*6*/ {{3, 4}, 30.0}, /*7*/ {{3, 5}, 80.0},
            /*8*/ {{4, 5}, 40.0}
        },
        mg_graph::GraphType::kDirectedGraph
    );
    std::vector<shortest_paths::Path<>> expected_paths = {
        {{0, 1, 3, 4, 5}, {0, 3, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        {{0, 2, 3, 4, 5}, {1, 4, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        {{0, 3, 4, 5}, {2, 6, 8}, {0.0, 100.0, 130.0, 170.0}, 170.0}
    };

    // Find the top 3 best paths
    auto paths = shortest_paths::KShortestPaths(*G, 0, 5, 3, check_abort_noop);

    // Could technically just compare the vectors, but this gives more useful output.
    ASSERT_EQ(paths.size(), expected_paths.size());
    ASSERT_TRUE(CheckPaths(paths, expected_paths));
}

/*
 *                ┌───┐                          ┌───┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └─▲─┘         │                └▲─┬┘
 *   │              │          40                 │ │
 *   ├──70──────────┘           │                 │ 40
 *   │                          │                 │ │
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     ├──80─────►│ 5 │
 *   │                          │     │          └─▲─┘
 *   │            ┌───┐         │     │            │
 *   └──50───────►│ 1 ├─────────┘     └──90────────┘
 *                └───┘
 */
TEST(ShortestPaths, YensSmallAcyclicGraphParallelEdges) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, 30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{0, 2}, 70.0},
            /*10*/ {{3, 5}, 90.0}
        },
        mg_graph::GraphType::kDirectedGraph
    );

    std::vector<shortest_paths::Path<>> expected_paths = {
        /* 0*/ {{0, 1, 3, 4, 5}, {0, 3, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        /* 1*/ {{0, 2, 3, 4, 5}, {1, 4, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        /* 2*/ {{0, 1, 3, 5}, {0, 3, 7}, {0.0, 50.0, 90.0, 170.0}, 170.0},
        /* 3*/ {{0, 3, 4, 5}, {2, 6, 8}, {0.0, 100.0, 130.0, 170.0}, 170.0},
        /* 4*/ {{0, 2, 4, 5}, {1, 5, 8}, {0.0, 50.0, 130.0, 170.0}, 170.0},
        /* 5*/ {{0, 2, 3, 5}, {1, 4, 7}, {0.0, 50.0, 90.0, 170.0}, 170.0},
        /* 6*/ {{0, 3, 5}, {2, 7}, {0.0, 100.0, 180.0}, 180.0},
        // New paths added by parallel edges
        /* 7*/ {{0, 2, 3, 5}, {1, 4, 10}, {0.0, 50.0, 90.0, 180.0}, 180.0},
        /* 8*/ {{0, 1, 3, 5}, {0, 3, 10}, {0.0, 50.0, 90.0, 180.0}, 180.0},
        /* 9*/ {{0, 2, 3, 4, 5}, {9, 4, 6, 8}, {0.0, 70.0, 110.0, 140.0, 180.0}, 180.0},
        /*10*/ {{0, 3, 5}, {2, 10}, {0.0, 100.0, 190.0}, 190.0},
        /*11*/ {{0, 2, 3, 5}, {9, 4, 7}, {0.0, 70.0, 110.0, 190.0}, 190.0},
        /*12*/ {{0, 2, 4, 5}, {9, 5, 8}, {0.0, 70.0, 150.0, 190.0}, 190.0},
        /*13*/ {{0, 2, 3, 5}, {9, 4, 10}, {0.0, 70.0, 110.0, 200.0}, 200.0}
    };

    // Set K really high so we find all the paths
    auto paths = shortest_paths::KShortestPaths(*G, 0, 5, 100, check_abort_noop);
    ASSERT_EQ(paths.size(), expected_paths.size());
    ASSERT_TRUE(CheckPaths(paths, expected_paths));
}

/*
 *                  ┌─────────────10───────────────┐
 *                ┌─▼─┐                          ┌─┴─┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, YensSmallGraphWithCycle) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, 30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{4, 2}, 10.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );

    // Should add no new paths over version without the cycle
    std::vector<shortest_paths::Path<>> expected_paths = {
        /* 0*/ {{0, 1, 3, 4, 5}, {0, 3, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        /* 1*/ {{0, 2, 3, 4, 5}, {1, 4, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        /* 2*/ {{0, 1, 3, 5}, {0, 3, 7}, {0.0, 50.0, 90.0, 170.0}, 170.0},
        /* 3*/ {{0, 3, 4, 5}, {2, 6, 8}, {0.0, 100.0, 130.0, 170.0}, 170.0},
        /* 4*/ {{0, 2, 4, 5}, {1, 5, 8}, {0.0, 50.0, 130.0, 170.0}, 170.0},
        /* 5*/ {{0, 2, 3, 5}, {1, 4, 7}, {0.0, 50.0, 90.0, 170.0}, 170.0},
        /* 6*/ {{0, 3, 5}, {2, 7}, {0.0, 100.0, 180.0}, 180.0},
    };

    // Set K really high so we find all the paths
    auto paths = shortest_paths::KShortestPaths(*G, 0, 5, 100, check_abort_noop);
    ASSERT_EQ(paths.size(), expected_paths.size());
    ASSERT_TRUE(CheckPaths(paths, expected_paths));
}

// Subset of a real graph
TEST(ShortestPaths, YensComplexGraph) {
    auto G = mg_generate::BuildWeightedGraph(
        21,
        {
            {{0, 1}, 44.35767007519712}, {{0, 1}, 44.3614195558365}, {{0, 2}, 36.96636525857858}, {{2, 1}, 24.169029158536663},
            {{2, 1}, 24.170524221282097}, {{2, 1}, 24.171429038544662}, {{0, 2}, 64.55024698269057}, {{0, 2}, 64.55141007312834},
            {{0, 2}, 64.55167609024404}, {{2, 1}, 51.75544726946053}, {{0, 2}, 64.55410993813733}, {{0, 3}, 63.242882527792375},
            {{3, 1}, 25.49204070545458}, {{0, 4}, 38.14445427698793}, {{4, 2}, 42.85208047898673}, {{0, 5}, 44.19882394401838},
            {{5, 1}, 44.5211018393918}, {{0, 6}, 42.59554957482579}, {{6, 1}, 46.11372734572857}, {{0, 7}, 58.35094417766534},
            {{7, 1}, 30.045636085948253}, {{0, 8}, 60.546326163193726}, {{8, 1}, 28.569578964946295}, {{0, 9}, 56.13689119148638},
            {{9, 1}, 32.57877915263091}, {{0, 10}, 46.08582581107471}, {{10, 1}, 42.661210705676844}, {{0, 11}, 49.90815863607187},
            {{11, 1}, 38.82112192850234}, {{0, 12}, 42.894284708141846}, {{12, 1}, 45.838855839588916}, {{6, 2}, 38.171378253852374},
            {{7, 2}, 22.666921251998378}, {{7, 2}, 22.69960262863917}, {{4, 2}, 42.97544625155755}, {{3, 2}, 17.904590234969124},
            {{11, 2}, 31.41363979609325}, {{5, 2}, 37.12298668941268}, {{10, 2}, 35.24417692059526}, {{10, 13}, 42.67973648278156},
            {{13, 1}, 16.730401083744233}, {{9, 2}, 25.19473940461525}, {{9, 2}, 25.202077506979744}, {{12, 2}, 38.45161504113953},
            {{12, 2}, 38.45234973891122}, {{8, 2}, 21.18158470954177}, {{0, 12}, 70.46381445783558}, {{0, 14}, 37.340202958769225},
            {{14, 1}, 51.38716813099613}, {{0, 13}, 44.3644818161463}, {{0, 13}, 71.99243802792876}, {{0, 15}, 48.947627139594886},
            {{15, 1}, 12.181628518565695}, {{0, 16}, 34.194599107966695}, {{16, 1}, 54.51447429966375}, {{0, 17}, 41.534585643389576},
            {{17, 1}, 47.00149287709005}, {{13, 1}, 44.36072669472187}, {{0, 18}, 42.06374971722343}, {{18, 1}, 52.02270384142229},
            {{0, 19}, 37.75827077384188}, {{19, 1}, 43.4637943430707}, {{0, 20}, 39.26248104943602}, {{20, 1}, 10.99298794961728},
            {{20, 2}, 31.27024623470097}, {{19, 2}, 35.83979778672537}, {{17, 2}, 35.790273413104174}, {{18, 2}, 38.830069447697404},
            {{13, 2}, 36.920091616035066}, {{13, 2}, 36.922191511042065}, {{14, 2}, 43.96983913453139}

        },
        mg_graph::GraphType::kDirectedGraph
    );
    const uint64_t K = 100;

    // Should be more than 100 paths
    auto paths = shortest_paths::KShortestPaths(*G, 0, 1, K, check_abort_noop);
    ASSERT_EQ(paths.size(), K);
    double prev_weight = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < paths.size(); i++) {
        const auto& path = paths[i];
        ASSERT_GE(path.total_cost, prev_weight) << "Path " << i << " had weight less than prevous path";
        prev_weight = path.total_cost;
    }

    // Get paths in single-threaded mode, should be the same
    auto paths_seq = shortest_paths::KShortestPaths(*G, 0, 1, K, check_abort_noop, 1);
    ASSERT_EQ(paths_seq.size(), paths.size());
    for (size_t i = 0; i < paths_seq.size(); i++) {
        const auto& path = paths_seq[i];
        ASSERT_EQ(path.edges, paths[i].edges) << "Path " << i << " has different edges than in parallel result";
    }
}

TEST(ShortestPaths, BellmanFordEmptyGraph) {
    auto G = mg_generate::BuildGraph(0, {});

    auto path = shortest_paths::BellmanFord<uint64_t>(*G, 0, 0, {}, {}, check_abort_noop);

    ASSERT_TRUE(path.empty());
}

TEST(ShortestPaths, BellmanFordSingleNode) {
    auto G = mg_generate::BuildGraph(1, {});

    auto path = shortest_paths::BellmanFord<uint64_t>(*G, 0, 0, {}, {}, check_abort_noop);

    ASSERT_TRUE(path.empty());
}

TEST(ShortestPaths, BellmanFordDisconnectedNodes) {
    auto G = mg_generate::BuildGraph(10, {});

    auto path = shortest_paths::BellmanFord<uint64_t>(*G, 0, 9, {}, {}, check_abort_noop);

    ASSERT_TRUE(path.empty());
}

/*
 *                  ┌─────────────10───────────────┐
 *                ┌─▼─┐                          ┌─┴─┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, BellmanFordSmallGraphWithCycle) {
    size_t NUM_VERTEX = 6;
    uint64_t SOURCE = 0;
    auto G = mg_generate::BuildWeightedGraph(
        NUM_VERTEX,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, 30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{4, 2}, 10.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );
    // Expected paths from 0 -> i
    std::vector<shortest_paths::Path<>> expected_paths{
        /*0->0*/ {{0}, {}, {0.0}, 0.0},
        /*0->1*/ {{0, 1}, {0}, {0.0, 50.0}, 50.0},
        /*0->2*/ {{0, 2}, {1}, {0.0, 50.0}, 50.0},
        /*0->3*/ {{0, 1, 3}, {0, 3}, {0.0, 50.0, 90.0}, 90.0},
        /*0->4*/ {{0, 1, 3, 4}, {0, 3, 6}, {0.0, 50.0, 90.0, 120.0}, 120.0},
        /*0->5*/ {{0, 1, 3, 4, 5}, {0, 3, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0}
    };

    shortest_paths::BellmanFordPathfinder<> pathfinder(*G, SOURCE);
    ASSERT_FALSE(pathfinder.has_negative_cycle());
    for (uint64_t target = 0; target  < NUM_VERTEX; target++) {
        ASSERT_TRUE(pathfinder.has_path_to(target));
        auto path = pathfinder.path_to(target);
        ASSERT_EQ(path, expected_paths[target]);
    }
}

/*
 *                  ┌─────────── -100 ─────────────┐
 *                ┌─▼─┐                          ┌─┴─┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬─-30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, BellmanFordSmallGraphWithNegativeCycle) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, -30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{4, 2}, -100.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );
    shortest_paths::Path<> expected_cycle{
        {2, 4, 2},
        {5, 9},
        {0.0, 80.0, -20.0},
        -20.0,
    };

    shortest_paths::BellmanFordPathfinder<> pathfinder(*G, 0);
    ASSERT_TRUE(pathfinder.has_negative_cycle());
    auto cycle = pathfinder.negative_cycle().value();
    ASSERT_EQ(cycle.nodes, expected_cycle.nodes);
    ASSERT_EQ(cycle.edges, expected_cycle.edges);
    ASSERT_EQ(cycle.costs, expected_cycle.costs);
    ASSERT_EQ(cycle.total_cost, expected_cycle.total_cost);
}

TEST(ShortestPaths, IterativeBellmanFordEmptyGraph) {
    auto G = mg_generate::BuildGraph(0, {});

    auto path = shortest_paths::IterativeBellmanFord<uint64_t>(*G, 0, 0, {}, {}, check_abort_noop);

    ASSERT_TRUE(path.empty());
}

TEST(ShortestPaths, IterativeBellmanFordSingleNode) {
    auto G = mg_generate::BuildGraph(1, {});

    auto path = shortest_paths::IterativeBellmanFord<uint64_t>(*G, 0, 0, {}, {}, check_abort_noop);

    ASSERT_TRUE(path.empty());
}

TEST(ShortestPaths, IterativeBellmanFordDisconnectedNodes) {
    auto G = mg_generate::BuildGraph(10, {});

    auto path = shortest_paths::IterativeBellmanFord<uint64_t>(*G, 0, 9, {}, {}, check_abort_noop);

    ASSERT_TRUE(path.empty());
}

/*
 *                  ┌─────────────10───────────────┐
 *                ┌─▼─┐                          ┌─┴─┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, IterativeBellmanFordSmallGraphWithCycle) {
    size_t NUM_VERTEX = 6;
    auto G = mg_generate::BuildWeightedGraph(
        NUM_VERTEX,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, 30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{4, 2}, 10.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );
    shortest_paths::Path<> expected_path = {{0, 1, 3, 4, 5}, {0, 3, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0};

    auto path = shortest_paths::IterativeBellmanFord<uint64_t>(*G, 0, 5, {}, {}, check_abort_noop);
    ASSERT_EQ(path, expected_path);
}

/*
 *                  ┌─────────── -100 ─────────────┐
 *                ┌─▼─┐                          ┌─┴─┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬─-30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, IterativeBellmanFordSmallGraphWithNegativeCycle) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, -30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{4, 2}, -100.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );
    std::vector<double> edge_scores;
    for (const auto& edge : G->Edges()) {
        edge_scores.push_back(G->GetWeight(edge.id));
    }

    // Should end up with the same path regardless of order.
    // For ascending, it should remove the -100 edge, which breaks the cycle.
    // For descending, it should remove both positive edges leading out of node 2, which makes
    // the -100 edge unusable, which also breaks the cycle.
    shortest_paths::Path<> expected_path{
        {0, 1, 3, 4, 5},
        {0, 3, 6, 8},
        {0.0, 50.0, 90.0, 60.0, 100.0},
        100.0
    };

    shortest_paths::IterativeBellmanFordPathfinder<uint64_t> pathfinder(edge_scores, true);
    auto path = pathfinder.search(*G, 0, 5, {}, {}, check_abort_noop);
    ASSERT_EQ(path, expected_path);
    ASSERT_EQ(pathfinder.num_edges_removed(), 1ULL);

    pathfinder.cull_ascending(false);
    path = pathfinder.search(*G, 0, 5, {}, {}, check_abort_noop);
    ASSERT_EQ(path, expected_path);
    ASSERT_EQ(pathfinder.num_edges_removed(), 2ULL);
}

TEST(ShortestPaths, JohnsonsEmptyGraph) {
    auto G = mg_generate::BuildGraph(0, {});

    shortest_paths::JohnsonsPathfinder<> pathfinder;
    pathfinder.search_all(*G, check_abort_noop, 1);

    ASSERT_EQ(pathfinder.num_reachable_nodes_from(0), 0UL);
}

TEST(ShortestPaths, JohnsonsSingleNode) {
    auto G = mg_generate::BuildGraph(1, {});

    shortest_paths::JohnsonsPathfinder<> pathfinder;
    pathfinder.search_all(*G, check_abort_noop, 1);

    ASSERT_EQ(pathfinder.num_reachable_nodes_from(0), 1UL);
}

TEST(ShortestPaths, JohnsonsDisconnectedNodes) {
    auto G = mg_generate::BuildGraph(10, {});

    shortest_paths::JohnsonsPathfinder<> pathfinder;
    pathfinder.search_all(*G, check_abort_noop, 1);

    ASSERT_EQ(pathfinder.num_reachable_nodes_from(0), 1UL);

    std::vector<std::uint64_t> sources = {0};
    pathfinder.search_some(*G, sources, check_abort_noop, 1);

    ASSERT_EQ(pathfinder.num_reachable_nodes_from(0), 1UL);
    ASSERT_FALSE(pathfinder.has_pathfinder(1));

    // Test some other container types and thread numbers to make sure all methods are instantiated.
    std::unordered_set<uint64_t> sources_uset = {0};
    std::set<uint64_t> sources_set = {0};
    std::list<uint64_t> sources_list = {0};
    std::forward_list<uint64_t> sources_flist = {0};

    pathfinder.search_some(*G, sources_uset, check_abort_noop, 2);
    ASSERT_EQ(pathfinder.num_reachable_nodes_from(0), 1UL);
    ASSERT_FALSE(pathfinder.has_pathfinder(1));

    std::vector<double> scores(10, 0.0);
    pathfinder.search_some_remove_cycles(*G, sources_set, scores, true, check_abort_noop, 1);
    ASSERT_EQ(pathfinder.num_reachable_nodes_from(0), 1UL);
    ASSERT_FALSE(pathfinder.has_pathfinder(1));

    pathfinder.search_some_remove_cycles(*G, sources_list, scores, true, check_abort_noop, 2);
    ASSERT_EQ(pathfinder.num_reachable_nodes_from(0), 1UL);
    ASSERT_FALSE(pathfinder.has_pathfinder(1));

    pathfinder.search_some(*G, sources_flist, check_abort_noop, 1);
    ASSERT_EQ(pathfinder.num_reachable_nodes_from(0), 1UL);
    ASSERT_FALSE(pathfinder.has_pathfinder(1));
}

/*
 *                  ┌────────────-10───────────────┐
 *                ┌─▼─┐                          ┌─┴─┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, JohnsonsSmallGraphWithNegativeEdge) {
    size_t NUM_VERTEX = 6;
    auto G = mg_generate::BuildWeightedGraph(
        NUM_VERTEX,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, 30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{4, 2}, -10.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );
    std::vector<double> expected_node_weights = {0.0, 0.0, -10.0, 0.0, 0.0, 0.0};
    std::vector<double> expected_edge_weights = {
        50.0, 60.0, 100.0, 40.0, 30.0, 70.0, 30.0, 80.0, 40.0, 0.0
    };

    shortest_paths::JohnsonsPathfinder<> pathfinder;
    pathfinder.search_all(*G, check_abort_noop, 0);

    ASSERT_FALSE(pathfinder.has_negative_cycle());
    ASSERT_EQ(pathfinder.node_weights(), expected_node_weights);
    ASSERT_EQ(pathfinder.edge_weights(), expected_edge_weights);

    for (size_t source_id = 0; source_id < NUM_VERTEX; source_id++) {
        ASSERT_TRUE(pathfinder.has_pathfinder(source_id));
    }

    ASSERT_TRUE(pathfinder.has_path(0, 5));
    auto path_0_5 = pathfinder.get_path(0, 5);
    ASSERT_FLOAT_EQ(path_0_5.total_cost, 160.0);

    auto path_2_3 = pathfinder.get_path(2, 3);
    ASSERT_FLOAT_EQ(path_2_3.total_cost, 30.0);
}

/*
 *                  ┌─────────── -100 ─────────────┐
 *                ┌─▼─┐                          ┌─┴─┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬─-30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, JohnsonsSmallGraphWithNegativeCycle) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, -30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{4, 2}, -100.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );
    std::vector<double> edge_scores;
    for (const auto& edge : G->Edges()) {
        edge_scores.push_back(G->GetWeight(edge.id));
    }
    shortest_paths::Path<> expected_cycle{
        {2,3,  4, 2},
        {4, 6, 9},
        {0.0, 40.0, 10.0, -90.0},
        -90.0,
    };

    shortest_paths::JohnsonsPathfinder<> pathfinder;
    pathfinder.search_all(*G, check_abort_noop, 0);

    ASSERT_TRUE(pathfinder.has_negative_cycle());
    ASSERT_EQ(pathfinder.negative_cycle().value(), expected_cycle);

    std::vector<double> expected_node_weights = {0.0, 0.0, 0.0, 0.0, -30.0, 0.0};
    std::vector<double> expected_edge_weights = {
        50.0, 50.0, 100.0, 40.0, 40.0, 110.0, 0.0, 80.0, 10.0, -130.0
    };
    std::unordered_set<uint64_t> expedted_removed_edges = {9};

    pathfinder.search_all_remove_cycles(*G, edge_scores, true, check_abort_noop, 0);
    ASSERT_FALSE(pathfinder.has_negative_cycle());
    ASSERT_EQ(pathfinder.num_edges_removed(), 1UL);
    ASSERT_EQ(pathfinder.removed_edges(), expedted_removed_edges);
    ASSERT_EQ(pathfinder.node_weights(), expected_node_weights);
    ASSERT_EQ(pathfinder.edge_weights(), expected_edge_weights);

    shortest_paths::Path<> expected_path_0_5 = {
        {0, 1, 3, 4, 5},
        {0, 3, 6, 8},
        {0.0, 50.0, 90.0, 90.0, 100.0},
        100.0
    };

    auto path_0_5 = pathfinder.get_path(0, 5);
    ASSERT_EQ(path_0_5.nodes, expected_path_0_5.nodes);
    ASSERT_EQ(path_0_5.edges, expected_path_0_5.edges);
    ASSERT_EQ(path_0_5.costs, expected_path_0_5.costs);
    ASSERT_FLOAT_EQ(path_0_5.total_cost, expected_path_0_5.total_cost);
}

/*
 *                  ┌─────────── -100 ─────────────┐
 *                ┌─▼─┐                          ┌─┴─┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬─-30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, JohnsonsKShortestSmallNegCycle) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /* 0*/ {{0, 1}, 50.0}, /* 1*/ {{0, 2}, 50.0}, /* 2*/ {{0, 3}, 100.0},
            /* 3*/ {{1, 3}, 40.0},
            /* 4*/ {{2, 3}, 40.0}, /* 5*/ {{2, 4}, 80.0},
            /* 6*/ {{3, 4}, -30.0}, /* 7*/ {{3, 5}, 80.0},
            /* 8*/ {{4, 5}, 40.0},
            // new edges
            /* 9*/ {{4, 2}, -100.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );
    std::vector<double> edge_scores;
    for (const auto& edge : G->Edges()) {
        edge_scores.push_back(G->GetWeight(edge.id));
    }

    // Expected paths once the cycle is eliminated
    std::vector<shortest_paths::Path<>> expected_paths = {
        {{0, 1, 3, 4, 5}, {0, 3, 6, 8}, {0.0, 50.0, 90.0, 90.0, 100.0}, 100.0},
        {{0, 2, 3, 4, 5}, {1, 4, 6, 8}, {0.0, 50.0, 90.0, 90.0, 100.0}, 100.0},
        {{0, 3, 4, 5}, {2, 6, 8}, {0.0, 100.0, 100.0, 110.0}, 110.0}
    };

    shortest_paths::JohnsonsPathfinder<> pathfinder;
    auto paths = pathfinder.k_shortest_paths(*G, 0, 5, 3, check_abort_noop);
    ASSERT_TRUE(paths.empty());
    ASSERT_TRUE(pathfinder.has_negative_cycle());

    // Find the top 3 best paths
    paths = pathfinder.k_shortest_paths_remove_cycles(*G, 0, 5, 3, edge_scores, true, check_abort_noop);

    // Could technically just compare the vectors, but this gives more useful output.
    ASSERT_EQ(paths.size(), expected_paths.size());
    ASSERT_TRUE(CheckPaths(paths, expected_paths));
}

TEST(ShortestPaths, DisjointEmptyGraph) {
    auto G = mg_generate::BuildGraph(0, {});

    auto paths = shortest_paths::DisjointKShortestPaths(*G, 0UL, 0UL, 0UL, check_abort_noop);

    ASSERT_TRUE(paths.empty());
}

/*
 *                ┌───┐                          ┌───┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, DisjointSmallAcyclicGraph) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /*0*/ {{0, 1}, 50.0}, /*1*/ {{0, 2}, 50.0}, /*2*/ {{0, 3}, 100.0},
            /*3*/ {{1, 3}, 40.0},
            /*4*/ {{2, 3}, 40.0}, /*5*/ {{2, 4}, 80.0},
            /*6*/ {{3, 4}, 30.0}, /*7*/ {{3, 5}, 80.0},
            /*8*/ {{4, 5}, 40.0}
        },
        mg_graph::GraphType::kDirectedGraph
    );
    std::vector<shortest_paths::Path<>> expected_paths = {
        {{0, 1, 3, 4, 5}, {0, 3, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        {{0, 2, 3, 5}, {1, 4, 7}, {0.0, 50.0, 90.0, 170.0}, 170.0},
    };

    // Find all disjoint shortest paths
    auto paths = shortest_paths::DisjointKShortestPaths<>(*G, 0UL, 5UL, 0UL, check_abort_noop);

    // Could technically just compare the vectors, but this gives more useful output.
    ASSERT_EQ(paths.size(), expected_paths.size());
    ASSERT_TRUE(CheckPaths(paths, expected_paths));
}

/*
 *                ┌───┐                          ┌───┐
 *   ┌──50───────►│ 2 ├─────────┬────────80─────►│ 4 │
 *   │            └───┘         │                └▲─┬┘
 *   │                         40                 │ │
 *   │                          │                 │ 40
 * ┌─┴─┐                       ┌▼──┐              │ │
 * │ 0 ├────────100───────────►│ 3 ├──┬──30───────┘ │
 * └─┬─┘                       └▲──┘  │             │
 *   │                          │     │          ┌──▼┐
 *   │                         40     └──80─────►│ 5 │
 *   │                          │                └───┘
 *   │            ┌───┐         │
 *   └──50───────►│ 1 ├─────────┘
 *                └───┘
 */
TEST(ShortestPaths, PartialDisjointSmallAcyclicGraph) {
    auto G = mg_generate::BuildWeightedGraph(
        6,
        {
            /*0*/ {{0, 1}, 50.0}, /*1*/ {{0, 2}, 50.0}, /*2*/ {{0, 3}, 100.0},
            /*3*/ {{1, 3}, 40.0},
            /*4*/ {{2, 3}, 40.0}, /*5*/ {{2, 4}, 80.0},
            /*6*/ {{3, 4}, 30.0}, /*7*/ {{3, 5}, 80.0},
            /*8*/ {{4, 5}, 40.0}
        },
        mg_graph::GraphType::kDirectedGraph
    );
    // Score edges by their ID number
    std::vector<double> scores = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<shortest_paths::Path<>> expected_paths = {
        {{0, 1, 3, 4, 5}, {0, 3, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        {{0, 2, 3, 4, 5}, {1, 4, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        {{0, 3, 4, 5}, {2, 6, 8}, {0.0, 100.0, 130.0, 170.0}, 170.0},
    }; 

    // Find all partial disjoint shortest paths, removing smallest-ID edge first
    auto paths = shortest_paths::PartialDisjointKShortestPaths<>(*G, 0UL, 5UL, 0UL, scores, 1UL, true, check_abort_noop);

    // Could technically just compare the vectors, but this gives more useful output.
    ASSERT_EQ(paths.size(), expected_paths.size());
    ASSERT_TRUE(CheckPaths(paths, expected_paths));
}

TEST(ShortestPaths, SSP_EmptyGraph) {
    auto G = mg_generate::BuildGraph(0, {});

    shortest_paths::SuccessiveShortestPathsPathfinder<uint64_t> pathfinder;
    std::vector<double> capacities;
    std::vector<double> factors;
    auto paths = pathfinder.search(
        *G, 0, 0, 0.0,
        capacities, factors,
        0.001, shortest_paths::FlowConversion::None, check_abort_noop
    );

    ASSERT_TRUE(paths.empty());
}

/* Edges marked as weight:capacity
 *
 *       ┌───2:30────┐
 *       │           │
 *       │           │
 * ┌───┐ │           │  ┌───┐               ┌───┐
 * │ 0 ├─┼───1:20────┼─►│ 1 ├─────1:110────►│ 2 │
 * └───┘ │           │  └───┘               └───┘
 *       │           │
 *       │           │
 *       └───3:60────┘
 */
TEST(ShortestPaths, SSP_SmallGraphNoConv) {
    auto G = mg_generate::BuildWeightedGraph(
        3,
        {
            /*0*/ {{0, 1}, 2.0},
            /*1*/ {{0, 1}, 1.0},
            /*2*/ {{0, 1}, 3.0},
            /*3*/ {{1, 2}, 1.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );

    using EdgesAndFlows = std::pair<std::vector<uint64_t>, std::vector<double>>;
    using EdgesAndFlowsVec = std::vector<EdgesAndFlows>;

    const EdgesAndFlowsVec expected_max_flow = {
        {/*edges*/ {1, 3}, /*flows*/ {20.0, 20.0, 20.0}},
        {/*edges*/ {0, 3}, /*flows*/ {30.0, 30.0, 30.0}},
        {/*edges*/ {2, 3}, /*flows*/ {60.0, 60.0, 60.0}},
    };

    const std::vector<double> capacities = {30.0, 20.0, 60.0, 110.0};
    const std::vector<double> factors; // Empty, not needed
    const double EPSILON = 1.0e-6;

    shortest_paths::SuccessiveShortestPathsPathfinder<uint64_t> pathfinder;
    // Search with flow_in = max flow for graph
    auto paths = pathfinder.search(
        *G, 0, 2, 110.0, capacities, factors,
        EPSILON, shortest_paths::FlowConversion::None, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_max_flow.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_max_flow[i].first);
        ASSERT_EQ(paths[i].second, expected_max_flow[i].second);
    }

    // Search with flow_in > max flow for graph, result should be identical to
    // previous.
    paths = pathfinder.search(
        *G, 0, 2, 100000.0, capacities, factors,
        EPSILON, shortest_paths::FlowConversion::None, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_max_flow.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_max_flow[i].first);
        ASSERT_EQ(paths[i].second, expected_max_flow[i].second);
    }

    // Search with flow_in = 100, should not fully utilize edge 2.
    const EdgesAndFlowsVec expected_flow_100 = {
        {/*edges*/ {1, 3}, /*flows*/ {20.0, 20.0, 20.0}},
        {/*edges*/ {0, 3}, /*flows*/ {30.0, 30.0, 30.0}},
        {/*edges*/ {2, 3}, /*flows*/ {50.0, 50.0, 50.0}},
    };

    paths = pathfinder.search(
        *G, 0, 2, 100.0, capacities, factors,
        EPSILON, shortest_paths::FlowConversion::None, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_flow_100.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_flow_100[i].first);
        ASSERT_EQ(paths[i].second, expected_flow_100[i].second);
    }

    // Search with flow_in = 50, should ignore edge 2 entirely
    const EdgesAndFlowsVec expected_flow_50 = {
        {/*edges*/ {1, 3}, /*flows*/ {20.0, 20.0, 20.0}},
        {/*edges*/ {0, 3}, /*flows*/ {30.0, 30.0, 30.0}},
    };

    paths = pathfinder.search(
        *G, 0, 2, 50.0, capacities, factors,
        EPSILON, shortest_paths::FlowConversion::None, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_flow_50.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_flow_50[i].first);
        ASSERT_EQ(paths[i].second, expected_flow_50[i].second);
    }
}

/* Edges marked as weight:capacity
 *
 *       ┌───2:30────┐
 *       │           │
 *       │           │
 * ┌───┐ │           │  ┌───┐               ┌───┐
 * │ 0 ├─┼───1:20────┼─►│ 1 ├─────1:110────►│ 2 │
 * └───┘ │           │  └───┘               └───┘
 *       │           │
 *       │           │
 *       └───3:60────┘
 * 
 * Tests for this one assume the weight acts as a conversion factor for flows
 * over that edge (# units out per unit in).
 * If a weight is 2, that means you get 2 units of output for each
 * unit of input.
 * 
 * Capacities are assumed to be in terms of the input to the edge.
 */
TEST(ShortestPaths, SSP_SmallGraphConv_TargetOverSource) {
    auto G = mg_generate::BuildWeightedGraph(
        3,
        {
            /*0*/ {{0, 1}, 2.0},
            /*1*/ {{0, 1}, 1.0},
            /*2*/ {{0, 1}, 3.0},
            /*3*/ {{1, 2}, 1.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );

    using EdgesAndFlows = std::pair<std::vector<uint64_t>, std::vector<double>>;
    using EdgesAndFlowsVec = std::vector<EdgesAndFlows>;

    const std::vector<double> capacities = {30.0, 20.0, 60.0, 110.0};
    const double EPSILON = 1.0e-6;
    const auto CONV_MODE = shortest_paths::FlowConversion::TargetOverSource;

    // Set factors to weights for convenience
    std::vector<double> factors(G->Edges().size());
    for (uint64_t i = 0; i < G->Edges().size(); i++) {
        factors[i] = G->GetWeight(i);
    }

    shortest_paths::SuccessiveShortestPathsPathfinder<uint64_t> pathfinder;

    // Search with flow_in = max flow for graph
    const EdgesAndFlowsVec expected_max_flow = {
        {/*edges*/ {1, 3}, /*flows*/ {20.0, 20.0, 20.0}},
        {/*edges*/ {0, 3}, /*flows*/ {30.0, 60.0, 60.0}},
        {/*edges*/ {2, 3}, /*flows*/ {10.0, 30.0, 30.0}},
    };

    double max_flow_in = 0.0;
    for (const auto& info : expected_max_flow) {
        max_flow_in += info.second[0];
    }
    auto paths = pathfinder.search(
        *G, 0, 2, max_flow_in, capacities, factors,
        EPSILON, CONV_MODE, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_max_flow.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_max_flow[i].first);
        ASSERT_EQ(paths[i].second, expected_max_flow[i].second);
    }

    // Search with flow_in > max flow for graph, result should be identical to
    // previous.
    paths = pathfinder.search(
        *G, 0, 2, max_flow_in * 2, capacities, factors,
        EPSILON, CONV_MODE, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_max_flow.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_max_flow[i].first);
        ASSERT_EQ(paths[i].second, expected_max_flow[i].second);
    }

    // Search with flow_in = 50, should ignore edge 2 entirely
    const EdgesAndFlowsVec expected_flow_2_paths = {
        {/*edges*/ {1, 3}, /*flows*/ {20.0, 20.0, 20.0}},
        {/*edges*/ {0, 3}, /*flows*/ {30.0, 60.0, 60.0}},
    };

    paths = pathfinder.search(
        *G, 0, 2, 50.0, capacities, factors,
        EPSILON, CONV_MODE, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_flow_2_paths.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_flow_2_paths[i].first);
        ASSERT_EQ(paths[i].second, expected_flow_2_paths[i].second);
    }
}

/* Edges marked as weight:capacity
 *
 *       ┌───2:30────┐
 *       │           │
 *       │           │
 * ┌───┐ │           │  ┌───┐               ┌───┐
 * │ 0 ├─┼───1:20────┼─►│ 1 ├─────2:110────►│ 2 │
 * └───┘ │           │  └───┘               └───┘
 *       │           │
 *       │           │
 *       └───3:60────┘
 * 
 * Tests for this one assume the weight acts as a conversion factor for flows
 * over that edge, (# units in per unit out).
 * If a weight is 2, that means you need 2 units of input for each
 * unit of output.
 * 
 * Capacities are assumed to be in terms of the input to the edge.
 */
TEST(ShortestPaths, SSP_SmallGraphConv_SourceOverTarget) {
    auto G = mg_generate::BuildWeightedGraph(
        3,
        {
            /*0*/ {{0, 1}, 2.0},
            /*1*/ {{0, 1}, 1.0},
            /*2*/ {{0, 1}, 3.0},
            /*3*/ {{1, 2}, 2.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );

    using EdgesAndFlows = std::pair<std::vector<uint64_t>, std::vector<double>>;
    using EdgesAndFlowsVec = std::vector<EdgesAndFlows>;

    const std::vector<double> capacities = {30.0, 20.0, 60.0, 110.0};
    const double EPSILON = 1.0e-6;
    const auto CONV_MODE = shortest_paths::FlowConversion::SourceOverTarget;

    // Set factors to weights for convenience
    std::vector<double> factors(G->Edges().size());
    for (uint64_t i = 0; i < G->Edges().size(); i++) {
        factors[i] = G->GetWeight(i);
    }

    shortest_paths::SuccessiveShortestPathsPathfinder<uint64_t> pathfinder;

    // Search with flow_in = max flow for graph
    const EdgesAndFlowsVec expected_max_flow = {
        {/*edges*/ {1, 3}, /*flows*/ {20.0, 20.0, 10.0}},
        {/*edges*/ {0, 3}, /*flows*/ {30.0, 15.0, 7.5}},
        {/*edges*/ {2, 3}, /*flows*/ {60.0, 20.0, 10.0}},
    };

    double max_flow_in = 0.0;
    for (const auto& info : expected_max_flow) {
        max_flow_in += info.second[0];
    }
    auto paths = pathfinder.search(
        *G, 0, 2, max_flow_in, capacities, factors,
        EPSILON, CONV_MODE, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_max_flow.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_max_flow[i].first);
        ASSERT_EQ(paths[i].second, expected_max_flow[i].second);
    }

    // Search with flow_in > max flow for graph, result should be identical to
    // previous.
    paths = pathfinder.search(
        *G, 0, 2, max_flow_in * 2, capacities, factors,
        EPSILON, CONV_MODE, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_max_flow.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_max_flow[i].first);
        ASSERT_EQ(paths[i].second, expected_max_flow[i].second);
    }

    // Search with flow_in = 95, should not fully utilize edge 2.
    const EdgesAndFlowsVec expected_flow_under_max = {
        {/*edges*/ {1, 3}, /*flows*/ {20.0, 20.0, 10.0}},
        {/*edges*/ {0, 3}, /*flows*/ {30.0, 15.0, 7.5}},
        {/*edges*/ {2, 3}, /*flows*/ {45.0, 15.0, 7.5}},
    };

    paths = pathfinder.search(
        *G, 0, 2, 95.0, capacities, factors,
        EPSILON, CONV_MODE, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_flow_under_max.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_flow_under_max[i].first);
        ASSERT_EQ(paths[i].second, expected_flow_under_max[i].second);
    }

    // Search with flow_in = 50, should ignore edge 2 entirely
    const EdgesAndFlowsVec expected_flow_2_paths = {
        {/*edges*/ {1, 3}, /*flows*/ {20.0, 20.0, 10.0}},
        {/*edges*/ {0, 3}, /*flows*/ {30.0, 15.0, 7.5}},
    };

    paths = pathfinder.search(
        *G, 0, 2, 50.0, capacities, factors,
        EPSILON, CONV_MODE, check_abort_noop
    );

    ASSERT_EQ(paths.size(), expected_flow_2_paths.size());
    for (size_t i = 0; i < paths.size(); i++) {
        ASSERT_EQ(paths[i].first.edges, expected_flow_2_paths[i].first);
        ASSERT_EQ(paths[i].second, expected_flow_2_paths[i].second);
    }
}

TEST(ShortestPaths, PolyEdge_Construct) {
    using shortest_paths::PolyEdge;

    PolyEdge default_edge;
    ASSERT_EQ(default_edge.id(), std::numeric_limits<uint64_t>::max());
    ASSERT_EQ(default_edge.coefficients().size(), 1);
    ASSERT_EQ(default_edge.coefficients()[0], 1);

    PolyEdge edge_1(1UL);
    ASSERT_EQ(edge_1.id(), 1UL);
    ASSERT_EQ(edge_1.coefficients().size(), 1);
    ASSERT_EQ(edge_1.coefficients()[0], 1);

    // y = x + 1
    Eigen::VectorXd x{{0.0, 1.0}};
    Eigen::VectorXd y{{1.0, 2.0}};
    PolyEdge edge_linear(2UL, x, y);
    ASSERT_EQ(edge_linear.id(), 2UL);
    ASSERT_EQ(edge_linear.coefficients().size(), 2);
    // coefficients should be (1.0, 1.0)
    ASSERT_DOUBLE_EQ(edge_linear.coefficients()(0), 1.0);
    ASSERT_DOUBLE_EQ(edge_linear.coefficients()(1), 1.0);

    std::vector<double> x_vec{0.0, 1.0};
    std::vector<double> y_vec{1.0, 2.0};
    PolyEdge edge_linear_std(3UL, x_vec, y_vec);
    ASSERT_EQ(edge_linear_std.id(), 3UL);
    ASSERT_EQ(edge_linear_std.coefficients().size(), 2);
    // coefficients should be (1.0, 1.0)
    ASSERT_DOUBLE_EQ(edge_linear_std.coefficients()(0), 1.0);
    ASSERT_DOUBLE_EQ(edge_linear_std.coefficients()(1), 1.0);
}

TEST(ShortestPaths, EdgeNetwork_Construct) {
    using PolyEdge = shortest_paths::PolyEdge<uint64_t>;
    using EdgeNetwork = shortest_paths::EdgeNetwork<uint64_t>;
    
    constexpr const uint64_t N_PATHS = 2UL;
    constexpr const uint64_t PATH_LEN = 2UL;

    std::vector<std::vector<PolyEdge>> poly_edge_list;
    for (uint64_t i = 0; i < N_PATHS; i++) {
        std::vector<PolyEdge> row;
        
        for (uint64_t j = 0; j < PATH_LEN; j++) {
            row.emplace_back((i * N_PATHS) + j);
        }

        poly_edge_list.push_back(row);
    }

    EdgeNetwork test_network(poly_edge_list);
    ASSERT_EQ(test_network.paths().size(), N_PATHS);
    for (uint64_t i = 0; i < N_PATHS; i++) {
        ASSERT_EQ(test_network.paths()[i].size(), PATH_LEN);
    }
}

/* Edges marked as weight:capacity
 *
 *       ┌───2:30────┐
 *       │           │
 *       │           │
 * ┌───┐ │           │  ┌───┐               ┌───┐
 * │ 0 ├─┼───1:20────┼─►│ 1 ├─────1:110────►│ 2 │
 * └───┘ │           │  └───┘               └───┘
 *       │           │
 *       │           │
 *       └───3:60────┘
 */
TEST(ShortestPaths, OptimizeFlow_Simple) {
    auto G = mg_generate::BuildWeightedGraph(
        3,
        {
            /*0*/ {{0, 1}, 2.0},
            /*1*/ {{0, 1}, 1.0},
            /*2*/ {{0, 1}, 3.0},
            /*3*/ {{1, 2}, 1.0},
        },
        mg_graph::GraphType::kDirectedGraph
    );
    const std::vector<double> capacities = {30.0, 20.0, 60.0, 110.0};

    std::vector<shortest_paths::Path<uint64_t>> paths {
        {{0, 1, 2}, {0, 3}, {2, 1}, {3}},
        {{0, 1, 2}, {1, 3}, {1, 1}, {2}},
    };

    auto result = shortest_paths::OptimizeFlows<>(*G, paths);

    ASSERT_TRUE(result.success);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
