# Add all module files related to Shortest Paths module

set(shortest_paths_src
    shortest_paths_module.cpp
    algorithm/dijkstra.cpp
    algorithm/yens.cpp
    util/conversions.cpp
    util/map.cpp
    util/options.cpp
)

add_query_module(shortest_paths 1 "${shortest_paths_src}")

target_link_libraries(shortest_paths PRIVATE mg_utility fmt::fmt highs::highs)
target_include_directories(shortest_paths
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} highs fmt
)

# Module tests
if (NOT MAGE_CUGRAPH_ENABLE)
    include(GoogleTest)
    set(shortest_paths_test_src
        shortest_paths_module_test.cpp
        algorithm/dijkstra.cpp
        algorithm/yens.cpp
    )

    add_executable(shortest_paths_test "${shortest_paths_test_src}")
    target_include_directories(shortest_paths_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(shortest_paths_test PRIVATE mg_utility mage_gtest fmt::fmt highs::highs)
    gtest_add_tests(TARGET shortest_paths_test)
endif()
