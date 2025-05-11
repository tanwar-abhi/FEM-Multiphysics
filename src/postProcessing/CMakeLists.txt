
set(postProcessCode
	postProcess.cpp		postProcess.hpp
	postProcessingInputs.cpp	postProcessingInputs.hpp
)

# Create a static library (*.so) of element code for linking with final executable
add_library(postProcessing STATIC ${postProcessCode})

# Add directory to the current static shared library
target_include_directories(postProcessing PUBLIC ${CMAKE_SOURCE_DIR}/Libraries/eigen3.3)
target_include_directories(postProcessing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(postProcessing PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_include_directories(postProcessing PUBLIC ${CMAKE_SOURCE_DIR}/src/element)
target_include_directories(postProcessing PUBLIC ${CMAKE_SOURCE_DIR}/src/preProcessing)
target_include_directories(postProcessing PUBLIC ${CMAKE_SOURCE_DIR}/src/readJsonInputs)

# target_link_libraries(postProcess element)