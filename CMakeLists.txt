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

SimplifyCT2.hpp
SimplifyCT2.cpp

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
constants.h
utl.h

ContourTreeData.hpp
ContourTreeData.cpp
)



add_library(reebgraph OBJECT ${ALL_SRCS})
set_target_properties(reebgraph PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_subdirectory(${CMAKE_SOURCE_DIR}/external/pybind11
        ${CMAKE_CURRENT_BINARY_DIR}/pybind11)        
    
pybind11_add_module(pyrg $<TARGET_OBJECTS:reebgraph> pymain.cpp)
target_link_libraries(pyrg PUBLIC pybind11::module ${PYTHON_LIBRARY})
set_target_properties(pyrg PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" )
target_link_libraries(pyrg ${ALL_LIBS})


option(PYRG_TESTS "Enable testing the PYRG module" OFF)
if(PYRG_TESTS)    
    add_subdirectory(tests)
endif(PYRG_TESTS)

if (WIN32)
    set(PYRG_MODULE_NAME "pyrg.pyd")
else (WIN32)
    set(PYRG_MODULE_NAME "pyrg.so")
endif (WIN32)

add_custom_command(TARGET pyrg POST_BUILD COMMAND ${CMAKE_COMMAND} -E
                copy $<TARGET_FILE:pyrg> ${CMAKE_CURRENT_SOURCE_DIR}/examples/${PYRG_MODULE_NAME})


