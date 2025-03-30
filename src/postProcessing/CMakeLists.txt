
set(postProcessCode
	postProcess.cpp		postProcess.hpp
	postProcessingInputs.cpp	postProcessingInputs.hpp
)

# Create a static library (*.so) of element code for linking with final executable
add_library(postProcessing STATIC ${postProcessCode})

# Add directory to the current static shared library
target_include_directories(postProcessing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
