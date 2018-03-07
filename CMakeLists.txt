project(LP_MP-Cut)
cmake_minimum_required(VERSION 2.8.12)

# C++11
add_compile_options(-std=c++17)

# compiler options
add_definitions(-march=native)

add_subdirectory(external/LP_MP)
#add_subdirectory(external/opengm)
add_subdirectory(external/MQLib)
add_subdirectory(external/andres_graph)

add_library(LP_MP-Cut INTERFACE)
target_include_directories(LP_MP-Cut INTERFACE .)
target_include_directories(LP_MP-Cut INTERFACE include)
target_include_directories(LP_MP-Cut INTERFACE external/opengm/include)
target_include_directories(LP_MP-Cut INTERFACE external/andres_graph/include)
target_include_directories(LP_MP-Cut INTERFACE external/PEGTL)
target_link_libraries(LP_MP-Cut INTERFACE LP_MP)


# for opengm
find_package(HDF5 1.8.15 REQUIRED)
add_library(HDF5 INTERFACE)
target_include_directories(HDF5 INTERFACE ${HDF5_INCLUDE_DIR})
target_link_libraries(HDF5 INTERFACE ${HDF5_LIBRARIES})
target_compile_definitions(HDF5 INTERFACE ${HDF5_DEFINITIONS})
message(STATUS ${HDF5_LIBRARIES})
message(STATUS ${HDF5_INCLUDE_DIR})

#enable_testing()
#add_subdirectory(test)
add_subdirectory(src)
