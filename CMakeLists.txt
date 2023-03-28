# Set minimum CMake version
cmake_minimum_required(VERSION 3.15)

# Set the project using CUDA as the source language
project(FEM1D LANGUAGES CUDA)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# Set default  build mode to Release
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

# Set the CUDA architecture to 61
set(CMAKE_CUDA_ARCHITECTURES 61)

# Set the source AND header files
file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu)
set(HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
file(GLOB_RECURSE HEADER_FILES ${HEADER_DIR}/*.cuh)
include_directories(${HEADER_DIR})

# Add the executable
add_executable(${PROJECT_NAME} ${SRC_FILES} ${HEADER_FILES})
set(${PROJECT_NAME}_INCLUDE_DIRS ${HEADER_DIR} CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(FILES ${HEADER_FILES} DESTINATION include)