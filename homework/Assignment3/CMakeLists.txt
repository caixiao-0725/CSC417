#include source code
file(GLOB H_FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
file(GLOB H_SHARED_FILES ${CMAKE_CURRENT_SOURCE_DIR}/shared/include/*.h)
file(GLOB CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB CPP_SHARED_FILES ${CMAKE_CURRENT_SOURCE_DIR}/shared/src/*.cpp)

#include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/shared/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${LIBIGL_INCLUDE_DIR})

#build executable 
add_executable(a3-finite-elements-3d main.cpp assignment_setup.h ${H_FILES} ${H_SHARED_FILES} ${CPP_FILES} ${CPP_SHARED_FILES})
target_link_libraries(a3-finite-elements-3d igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)