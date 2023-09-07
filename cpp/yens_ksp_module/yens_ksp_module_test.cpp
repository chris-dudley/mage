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
        {{0, 1, 3, 5}, {0, 3, 7}, {0.0, 50.0, 90.0, 170.0}, 170.0}
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

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}