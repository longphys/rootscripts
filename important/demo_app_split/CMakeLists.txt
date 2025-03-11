cmake_minimum_required(VERSION 3.10)
project(demo)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# Find the ROOT package
find_package(ROOT REQUIRED)
if(NOT ROOT_FOUND)
  message(FATAL_ERROR "ROOT is required but was not found. Please install ROOT and set up your environment properly.")
endif()

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})
# Include project headers
include_directories(${PROJECT_SOURCE_DIR}/include)

# Locating headers and sources
file(GLOB sources ${PROJECT_SOURCE_DIR}/source/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Copy all mac files
file(GLOB MACRO_FILES "*.mac")
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

# Create an executable
add_executable(demo demo.cc)

# Link against ROOT libraries
target_link_libraries(demo ${ROOT_LIBRARIES})