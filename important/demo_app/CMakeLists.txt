cmake_minimum_required(VERSION 3.10)
project(DemoRootApp)

# Find the ROOT package
find_package(ROOT REQUIRED)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

# Copy all mac files
file(GLOB MACRO_FILES "*.mac")
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

# Create an executable
add_executable(demo demo.C)

# Link against ROOT libraries
target_link_libraries(demo ${ROOT_LIBRARIES})

