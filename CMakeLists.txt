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
add_definitions(-DYGL_OPENGL=1)

find_package(OpenGL REQUIRED)
if(APPLE)
    include_directories(/usr/local/include)
    link_directories(/usr/local/lib)
    find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/local/lib)
endif(APPLE)
if(WIN32)
    include_directories(${CMAKE_SOURCE_DIR}/apps/w32/include)
    link_directories(${CMAKE_SOURCE_DIR}/apps/w32/lib-vc2015)
    find_library(GLEW_LIBRARIES NAMES glew32 PATHS ${CMAKE_SOURCE_DIR}/apps/w32/lib-vc2015)
    find_library(GLFW_LIBRARY NAMES glfw3dll PATHS ${CMAKE_SOURCE_DIR}/apps/w32/lib-vc2015)
endif(WIN32)
if(UNIX AND NOT APPLE)
    include_directories(/usr/include /usr/local/include)
    find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/lib /usr/local/lib64 /usr/lib64 /usr/local/lib /usr/lib/x86_64-linux-gnu)
    find_package(GLEW REQUIRED)
endif(UNIX AND NOT APPLE)
add_library(yocto_gl src/yocto/yocto_gl.h src/yocto/yocto_gl.cpp src/yocto/ext/stb_image.cpp src/yocto/ext/nanosvg.cpp src/yocto/ext/imgui/imgui.cpp src/yocto/ext/imgui/imgui_draw.cpp src/yocto/ext/imgui/imgui_impl_glfw_gl3.cpp src/yocto/ext/imgui/imgui_extra_fonts.cpp)
target_link_libraries(yocto_gl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GLEW_LIBRARIES} ${X11_LIBRARIES} X11 Xxf86vm Xrandr GL GLU Xi dl Xinerama rt GLEW m Xcursor)

if(UNIX AND NOT APPLE)
    find_package(Threads REQUIRED)
    target_link_libraries(yocto_gl Threads::Threads)
endif(UNIX AND NOT APPLE)

add_library(voro++ src/voro++/src/voro++.hh src/voro++/src/c_loops.hh src/voro++/src/cell.hh
        src/voro++/src/common.hh src/voro++/src/config.hh src/voro++/src/container.hh
        src/voro++/src/container_prd.hh src/voro++/src/rad_option.hh src/voro++/src/pre_container.hh
        src/voro++/src/unitcell.hh src/voro++/src/v_base.hh src/voro++/src/v_compute.hh
        src/voro++/src/wall.hh src/voro++/src/worklist.hh src/voro++/src/voro++.cc
        src/voro++/src/c_loops.cc src/voro++/src/cell.cc src/voro++/src/common.cc
        src/voro++/src/container.cc src/voro++/src/container_prd.cc src/voro++/src/pre_container.cc
        src/voro++/src/unitcell.cc src/voro++/src/v_base.cc src/voro++/src/v_compute.cc
        src/voro++/src/wall.cc)
set_target_properties(voro++ PROPERTIES LINKER_LANGUAGE CXX)

add_executable(tree-gen src/tree-gen.cpp)

target_link_libraries(tree-gen yocto_gl voro++)