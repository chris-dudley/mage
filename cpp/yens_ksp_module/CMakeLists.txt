# Add all module files related to Yen's algorithm module
set(yens_ksp_src
    yens_ksp_module.cpp
    algorithm/dijkstra.cpp
    algorithm/yens.cpp)

add_query_module(yens_ksp 1 "${yens_ksp_src}")

target_link_libraries(yens_ksp PRIVATE mg_utility)
target_include_directories(yens_ksp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Module tests
if (NOT MAGE_CUGRAPH_ENABLE)
    include(GoogleTest)
    set(yens_ksp_test_src
        yens_ksp_module_test.cpp
        algorithm/dijkstra.cpp
        algorithm/yens.cpp)
    
    add_executable(yens_ksp_test "${yens_ksp_test_src}")
    target_include_directories(yens_ksp_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(yens_ksp_test PRIVATE mg_utility mage_gtest)
    gtest_add_tests(TARGET yens_ksp_test)
endif()