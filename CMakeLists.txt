# Define minimum version of CMAKE required for build process
cmake_minimum_required(VERSION 3.20.0)

# Set Project name (executable will be of same name)
project(simulationMPS VERSION 0.8 LANGUAGES C CXX)

# Path where executable will be generated (Optional, otherise will be in build folder)
# set(EXECUTABLE_OUTPUT_PATH ../)

# Set c++ version and standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
message(STATUS "## C++17 standard selected")

# Compile flag for optimisation or debugger
# set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "-g")

# Pthread library for multithreading 
find_package(Threads)
if (Threads_FOUND)
    # Compilation options
    add_compile_options(-pthread)
    message(STATUS "## Pthread (Multithreading) library found")
endif()


# Find installation of CppJson parser library
find_package(jsoncpp REQUIRED)
if (jsoncpp_FOUND)
    message(STATUS "## cppJSON library found")
endif()


find_package(spdlog 1.15.3 QUIET)
if (spdlog_FOUND)
    message(STATUS "## spdlog v1.15.3 logger library found in local")

elseif(NOT spdlog_FOUND)
    include(FetchContent)

    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.15.3 # pin a release tag never a branch
        GIT_SHALLOW     TRUE
    )

    FetchContent_MakeAvailable(spdlog)
endif()


# Add subdirectories containing other source files
add_subdirectory(src)

# # Treat all warnings as errors
# add_compile_options(-Wall)
# message(STATUS "## Enable warrnings in code")
# add_compile_options( -Werror)
# message(STATUS "## Treat all warnings as errors")


# Create executable
add_executable(${PROJECT_NAME} main.cpp)


target_link_libraries(${PROJECT_NAME} JsonCpp::JsonCpp)
target_link_libraries(${PROJECT_NAME} jsonInputs)
target_link_libraries(${PROJECT_NAME} materialsProperty)
# target_link_libraries(${PROJECT_NAME} element)
# target_link_libraries(${PROJECT_NAME} postProcessing)
# target_link_libraries(${PROJECT_NAME} preProcessing)
# target_link_libraries(${PROJECT_NAME} srcFunction)

# Link your executable to spdlog : This automatically handles include directories and compiler flags
target_link_libraries(${PROJECT_NAME} spdlog::spdlog_header_only)