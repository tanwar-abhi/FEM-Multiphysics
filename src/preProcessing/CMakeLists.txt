
set(preProcessingCode
	readMesh.cpp	readMesh.hpp
)


# Create a static library (*.so) of element code for linking with final executable
add_library(preProcessing STATIC ${preProcessingCode})

# Add directory to the current static shared library
target_include_directories(preProcessing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})