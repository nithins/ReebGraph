#cmake_minimum_required(VERSION 3.9)
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project (reebgraph LANGUAGES CXX)


SET(ALL_SRCS
ContourTree.hpp
ContourTree.cpp

TriMesh.hpp
TriMesh.cpp

TopologicalFeatures.hpp
TopologicalFeatures.cpp

test.hpp

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



add_library(reebgraph SHARED ${ALL_SRCS})
target_compile_features(reebgraph PUBLIC cxx_std_11)
set_target_properties(reebgraph PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(reebgraph ${ALL_LIBS})