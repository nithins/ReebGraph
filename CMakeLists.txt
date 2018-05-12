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


if (WIN32)
    set(PYRG_MODULE_NAME "pyrg.pyd")
else (WIN32)
    set(PYRG_MODULE_NAME "pyrg.so")
endif (WIN32)


option(PYRG_TESTS "Enable testing the PYRG module" ON)

if(PYRG_TESTS)
    add_custom_command(TARGET pyrg POST_BUILD COMMAND ${CMAKE_COMMAND} -E 
                    copy $<TARGET_FILE:pyrg> ${CMAKE_CURRENT_SOURCE_DIR}/tests/${PYRG_MODULE_NAME})

    set(PYTEST_FILES test.py)


    add_custom_target(check COMMAND ${PYTHON_EXECUTABLE} -m pytest ${PYTEST_FILES}
                  DEPENDS pyrg
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests ${PYBIND11_USES_TERMINAL})




endif(PYRG_TESTS)



