cmake_minimum_required(VERSION 3.9 FATAL_ERROR)




#if (WIN32)
#    set(PYRG_MODULE_NAME "pyrg.pyd")
#else (WIN32)
#    set(PYRG_MODULE_NAME "pyrg.so")
#endif (WIN32)

# add_executable(testmain testMain.cpp test.hpp)



pybind11_add_module(pyrgtest $<TARGET_OBJECTS:reebgraph> pytest.cpp test.hpp)
target_link_libraries(pyrgtest PUBLIC pybind11::module ${PYTHON_LIBRARY})
set_target_properties(pyrgtest PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" )
set_target_properties(pyrgtest PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )


#add_custom_command(TARGET pyrg POST_BUILD COMMAND ${CMAKE_COMMAND} -E
#                copy $<TARGET_FILE:pyrg> ${CMAKE_CURRENT_SOURCE_DIR}/${PYRG_MODULE_NAME})

set(PYTEST_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/test.py
    )


add_custom_target(check COMMAND ${PYTHON_EXECUTABLE} -m pytest ${PYTEST_FILES}
              DEPENDS pyrg pyrgtest
              WORKING_DIRECTORY ${CMAKE_BINARY_DIR} ${PYBIND11_USES_TERMINAL})
