#cmake_minimum_required(VERSION 3.9)
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project (reebgraph LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 14)


SET(ALL_SRCS
src/ContourTree.hpp
src/ContourTree.cpp

src/TriMesh.hpp
src/TriMesh.cpp

src/TopologicalFeatures.hpp
src/TopologicalFeatures.cpp

#test.hpp

src/SimplifyCT.hpp
src/SimplifyCT.cpp

src/SimplifyCT2.hpp
src/SimplifyCT2.cpp

src/SimFunction.hpp
src/ScalarFunction.hpp

src/Persistence.hpp
src/Persistence.cpp

src/MergeTree.hpp
src/MergeTree.cpp

#main.cpp

src/HyperVolume.hpp
src/HyperVolume.cpp

src/Grid3D.hpp
src/Grid3D.cpp

src/DisjointSets.hpp
src/constants.h
src/utl.h

src/ContourTreeData.hpp
src/ContourTreeData.cpp
)



add_library(reebgraph OBJECT ${ALL_SRCS})
set_target_properties(reebgraph PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_subdirectory(${CMAKE_SOURCE_DIR}/external/pybind11
        ${CMAKE_CURRENT_BINARY_DIR}/pybind11)        
    
pybind11_add_module(pyrg $<TARGET_OBJECTS:reebgraph> src/pymain.cpp)
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

add_custom_command(TARGET pyrg POST_BUILD COMMAND ${CMAKE_COMMAND} -E
                copy $<TARGET_FILE:pyrg> ${CMAKE_CURRENT_SOURCE_DIR}/ReebAnatomy/${PYRG_MODULE_NAME})

