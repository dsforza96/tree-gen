cmake_minimum_required (VERSION 3.5)

project (tree-gen)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# compile flags
if(APPLE)
    add_definitions(-Wno-missing-braces)
endif(APPLE)
if(MSVC)
    add_definitions(/D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
    set(CMAKE_CXX_FLAGS "/EHsc")
endif(MSVC)
add_definitions(-DYGL_OPENGL=0)

add_library(yocto_gl src/yocto/yocto_gl.h src/yocto/yocto_gl.cpp src/yocto/ext/stb_image.cpp src/yocto/ext/nanosvg.cpp)

if(UNIX AND NOT APPLE)
    find_package(Threads REQUIRED)
    target_link_libraries(yocto_gl Threads::Threads)
endif(UNIX AND NOT APPLE)

add_library(voro++ src/voro++-0.4.6/src/voro++.hh src/voro++-0.4.6/src/c_loops.hh src/voro++-0.4.6/src/cell.hh
        src/voro++-0.4.6/src/common.hh src/voro++-0.4.6/src/config.hh src/voro++-0.4.6/src/container.hh
        src/voro++-0.4.6/src/container_prd.hh src/voro++-0.4.6/src/rad_option.hh src/voro++-0.4.6/src/pre_container.hh
        src/voro++-0.4.6/src/unitcell.hh src/voro++-0.4.6/src/v_base.hh src/voro++-0.4.6/src/v_compute.hh
        src/voro++-0.4.6/src/wall.hh src/voro++-0.4.6/src/worklist.hh src/voro++-0.4.6/src/voro++.cc
        src/voro++-0.4.6/src/c_loops.cc src/voro++-0.4.6/src/cell.cc src/voro++-0.4.6/src/common.cc
        src/voro++-0.4.6/src/container.cc src/voro++-0.4.6/src/container_prd.cc src/voro++-0.4.6/src/pre_container.cc
        src/voro++-0.4.6/src/unitcell.cc src/voro++-0.4.6/src/v_base.cc src/voro++-0.4.6/src/v_compute.cc
        src/voro++-0.4.6/src/wall.cc)
set_target_properties(voro++ PROPERTIES LINKER_LANGUAGE CXX)

add_executable(tree-gen src/tree-gen.cpp)

target_link_libraries(tree-gen yocto_gl voro++)