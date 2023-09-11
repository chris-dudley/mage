#include <map>

#include <gtest/gtest.h>
#include <mg_generate.hpp>
#include <mg_graph.hpp>
#include <mg_test_utils.hpp>

#include "algorithm/shortest_path.hpp"
#include "algorithm/yens.hpp"

void check_abort_noop() {}

bool CheckPaths(const std::vector<yens_alg::Path<>>& paths, const std::vector<yens_alg::Path<>>& expected) {
    std::multimap<double, size_t> cost_to_path_index;
    for (size_t i = 0; i < expected.size(); i++) {
        cost_to_path_index.emplace(expected[i].total_weight, i);
    }

    // Paths of equal length could be returned in any order and still be valid
    for (const auto& path : paths) {
        auto range = cost_to_path_index.equal_range(path.total_weight);
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

TEST(YensKSP, DijkstraEmptyGraph) {
    auto G = mg_generate::BuildGraph(0, {});

    auto path = yens_alg::Dijkstra(*G, 0, 0, {}, {});

    ASSERT_TRUE(path.empty());
}

TEST(YensKSP, DijkstraSingleNode) {
    auto G = mg_generate::BuildGraph(1, {});

    auto path = yens_alg::Dijkstra(*G, 0, 0, {}, {});

    ASSERT_TRUE(path.empty());
}

TEST(YensKSP, DijkstraDisconnectedNodes) {
    auto G = mg_generate::BuildGraph(10, {});

    auto path = yens_alg::Dijkstra(*G, 0, 9, {}, {});

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
TEST(YensKSP, DijkstraSmallAcyclicGraph) {
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
    yens_alg::Path<> expected_path{
        {0, 1, 3, 4, 5}, // nodes
        {0, 3, 6, 8},    // edges
        {0.0, 50.0, 90.0, 120.0, 160.0}, // weights
        160.0 // total weight
    };

    auto path = yens_alg::Dijkstra(*G, 0, 5, {}, {});
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
TEST(YensKSP, DijkstraSmallAcyclicGraphParallelEdges) {
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
    yens_alg::Path<> expected_path{
        {0, 1, 3, 4, 5}, // nodes
        {0, 3, 6, 8},    // edges
        {0.0, 50.0, 90.0, 120.0, 160.0}, // weights
        160.0 // total weight
    };

    auto path = yens_alg::Dijkstra(*G, 0, 5, {}, {});
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
TEST(YensKSP, DijkstraSmallGraphWithCycle) {
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
    yens_alg::Path<> expected_path{
        {0, 1, 3, 4, 5}, // nodes
        {0, 3, 6, 8},    // edges
        {0.0, 50.0, 90.0, 120.0, 160.0}, // weights
        160.0 // total weight
    };

    auto path = yens_alg::Dijkstra(*G, 0, 5, {}, {});
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
TEST(YensKSP, DijkstraSmallGraphWithNegativeCycle) {
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
    yens_alg::Path<> expected_path{
        {0, 1, 3, 4, 5}, // nodes
        {0, 3, 6, 8},    // edges
        {0.0, 50.0, 90.0, 60.0, 100.0}, // weights
        100.0 // total weight
    };

    auto path = yens_alg::Dijkstra(*G, 0, 5, {}, {});
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
TEST(YensKSP, DijkstraAllNegativeCycles) {
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
    auto path = yens_alg::Dijkstra(*G, 0, 5, {}, {});
    ASSERT_FALSE(path.empty());
    ASSERT_EQ(path.size(), 5);
}

TEST(YensKSP, DijkstraHugeCyclicGraph) {
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
    auto path = yens_alg::Dijkstra(*G, 0, NUM_VERTICIES-1, {}, {});
    ASSERT_EQ(path.empty(), false);
    ASSERT_EQ(path.size(), SHORTEST_PATH_EDGES);
    ASSERT_DOUBLE_EQ(path.total_weight, (SHORTEST_PATH_EDGES) * 1.0);
}

TEST(YensKSP, YensEmptyGraph) {
    auto G = mg_generate::BuildGraph(0, {});

    auto paths = yens_alg::KShortestPaths(*G, 0, 0, 1, yens_alg::Dijkstra, check_abort_noop);

    ASSERT_TRUE(paths.empty());
}

TEST(YensKSP, YensSingleNode) {
    auto G = mg_generate::BuildGraph(1, {});

    auto paths = yens_alg::KShortestPaths(*G, 0, 0, 1, yens_alg::Dijkstra, check_abort_noop);

    ASSERT_TRUE(paths.empty());
}

TEST(YensKSP, YensDisconnectedNodes) {
    auto G = mg_generate::BuildGraph(10, {});

    auto paths = yens_alg::KShortestPaths(*G, 0, 9, 1, yens_alg::Dijkstra, check_abort_noop);

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
TEST(YensKSP, YensSmallAcyclicGraph) {
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
    std::vector<yens_alg::Path<>> expected_paths = {
        {{0, 1, 3, 4, 5}, {0, 3, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        {{0, 2, 3, 4, 5}, {1, 4, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        {{0, 3, 4, 5}, {2, 6, 8}, {0.0, 100.0, 130.0, 170.0}, 170.0}
    };

    // Find the top 3 best paths
    auto paths = yens_alg::KShortestPaths(*G, 0, 5, 3, yens_alg::Dijkstra, check_abort_noop);

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
TEST(YensKSP, YensSmallAcyclicGraphParallelEdges) {
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

    std::vector<yens_alg::Path<>> expected_paths = {
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
    auto paths = yens_alg::KShortestPaths(*G, 0, 5, 100, yens_alg::Dijkstra, check_abort_noop);
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
TEST(YensKSP, YensSmallGraphWithCycle) {
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
    std::vector<yens_alg::Path<>> expected_paths = {
        /* 0*/ {{0, 1, 3, 4, 5}, {0, 3, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        /* 1*/ {{0, 2, 3, 4, 5}, {1, 4, 6, 8}, {0.0, 50.0, 90.0, 120.0, 160.0}, 160.0},
        /* 2*/ {{0, 1, 3, 5}, {0, 3, 7}, {0.0, 50.0, 90.0, 170.0}, 170.0},
        /* 3*/ {{0, 3, 4, 5}, {2, 6, 8}, {0.0, 100.0, 130.0, 170.0}, 170.0},
        /* 4*/ {{0, 2, 4, 5}, {1, 5, 8}, {0.0, 50.0, 130.0, 170.0}, 170.0},
        /* 5*/ {{0, 2, 3, 5}, {1, 4, 7}, {0.0, 50.0, 90.0, 170.0}, 170.0},
        /* 6*/ {{0, 3, 5}, {2, 7}, {0.0, 100.0, 180.0}, 180.0},
    };

    // Set K really high so we find all the paths
    auto paths = yens_alg::KShortestPaths(*G, 0, 5, 100, yens_alg::Dijkstra, check_abort_noop);
    ASSERT_EQ(paths.size(), expected_paths.size());
    ASSERT_TRUE(CheckPaths(paths, expected_paths));
}

// Subset of a real graph
TEST(YensKSP, YensComplexGraph) {
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

    // Should be more than 30 paths
    auto paths = yens_alg::KShortestPaths(*G, 0, 1, K, yens_alg::Dijkstra, check_abort_noop);
    ASSERT_EQ(paths.size(), K);
    double prev_weight = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < paths.size(); i++) {
        const auto& path = paths[i];
        ASSERT_GE(path.total_weight, prev_weight) << "Path " << i << " had weight less than prevous path";
        prev_weight = path.total_weight;
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}