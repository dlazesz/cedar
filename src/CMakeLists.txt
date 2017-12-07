include_directories(${PROJECT_SOURCE_DIR}/src/)
SET(HEADERS ${PROJECT_SOURCE_DIR}/src/cedar.h ${PROJECT_SOURCE_DIR}/src/cedarpp.h)
SET(EXECUTABLES ${PROJECT_BINARY_DIR}/src/cedar ${PROJECT_BINARY_DIR}/src/mkcedar)

add_executable(cedar ${HEADERS} cedar.cc)
add_executable(mkcedar ${HEADERS} mkcedar.cc)
add_executable(simple ${HEADERS} simple.cc)

INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(PROGRAMS ${EXECUTABLES} DESTINATION bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

# uninstall target from: https://cmake.org/Wiki/CMake_FAQ#Can_I_do_.22make_uninstall.22_with_CMake.3F
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
