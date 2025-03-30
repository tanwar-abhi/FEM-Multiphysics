
set(readInputsCode
	partInput.cpp		partInput.hpp
	solverInputs.cpp	solverInputs.hpp
)

# Create a static library (*.so) of element code for linking with final executable
add_library(inputsJson STATIC ${readInputsCode})

# Add directory to the current static shared library
target_include_directories(inputsJson PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})