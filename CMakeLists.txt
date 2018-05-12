#cmake_minimum_required(VERSION 3.9)
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project (reebgraph LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 14)


SET(ALL_SRCS
ContourTree.hpp
ContourTree.cpp

TriMesh.hpp
TriMesh.cpp

TopologicalFeatures.hpp
TopologicalFeatures.cpp

#test.hpp

SimplifyCT.hpp
SimplifyCT.cpp

SimFunction.hpp
ScalarFunction.hpp

Persistence.hpp
Persistence.cpp

MergeTree.hpp
MergeTree.cpp

#main.cpp

HyperVolume.hpp
HyperVolume.cpp

Grid3D.hpp
Grid3D.cpp

DisjointSets.hpp

ContourTreeData.hpp
ContourTreeData.cpp
)



#add_library(reebgraph SHARED ${ALL_SRCS})

add_subdirectory(${CMAKE_SOURCE_DIR}/external/pybind11
        ${CMAKE_CURRENT_BINARY_DIR}/pybind11)        
    
pybind11_add_module(pyrg ${ALL_SRCS} pymain.cpp)
target_link_libraries(pyrg PUBLIC pybind11::module ${PYTHON_LIBRARY})
set_target_properties(pyrg PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" )
target_compile_features(pyrg PUBLIC cxx_std_14)
set_target_properties(pyrg PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(pyrg ${ALL_LIBS})


option(PYRG_TESTS "Enable testing the PYRG module" ON)
if(PYRG_TESTS)    
    add_subdirectory(tests)
endif(PYRG_TESTS)



