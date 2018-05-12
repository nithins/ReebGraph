cmake_minimum_required(VERSION 3.9 FATAL_ERROR)




if (WIN32)
    set(PYRG_MODULE_NAME "pyrg.pyd")
else (WIN32)
    set(PYRG_MODULE_NAME "pyrg.so")
endif (WIN32)

add_executable(testmain testMain.cpp)


#add_custom_command(TARGET pyrg POST_BUILD COMMAND ${CMAKE_COMMAND} -E
#                copy $<TARGET_FILE:pyrg> ${CMAKE_CURRENT_SOURCE_DIR}/${PYRG_MODULE_NAME})

set(PYTEST_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/test.py
    )


add_custom_target(check COMMAND ${PYTHON_EXECUTABLE} -m pytest ${PYTEST_FILES}
              DEPENDS pyrg
              WORKING_DIRECTORY ${CMAKE_BINARY_DIR} ${PYBIND11_USES_TERMINAL})
