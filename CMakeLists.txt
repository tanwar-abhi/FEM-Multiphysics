# Define minimum version of CMAKE required for build process
cmake_minimum_required(VERSION 3.12.0)

# Set Project name (executable will be of same name)
project(simulationMPS VERSION 0.5 LANGUAGE C CXX)

# Path where executable will be generated (Optional, otherise will be in build folder)
# set(EXECUTABLE_OUTPUT_PATH ../)

# Set c++ version and standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
message(STATUS "## C++17 standard selected")

# Compile flag for optimisation or debugger
set(CMAKE_CXX_FLAGS "-O3")
# set(CMAKE_CXX_FLAGS "-g")

# Pthread library for multithreading 
find_package(Threads)
if (Thread_FOUND)
    # Compilation options
    add_compile_options(-pthread)
    message(STATUS "## Pthread (Multithreading) library found")
endif()


# Find installation of CppJson parser library
find_package(jsoncpp REQUIRED)
if (jsoncpp_FOUND)
    message(STATUS "## cppJSON library found")
endif()


# Add subdirectories containing other source files
add_subdirectory(${CMAKE_SOURCE_DIR}/src)

# # Treat all warnings as errors
# add_compile_options(-Wall)
# message(STATUS "## Enable warrnings in code")
# add_compile_options( -Werror)
# message(STATUS "## Treat all warnings as errors")



