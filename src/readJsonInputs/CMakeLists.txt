
set(readInputsCode
	partInputs.cpp		partInputs.hpp
	solverInputs.cpp	solverInputs.hpp
	boundaryInputs.cpp	boundaryInputs.hpp
)

# Create a static library (*.so) of element code for linking with final executable
add_library(inputsJson STATIC ${readInputsCode})

# Include directories for jsoncpp
# include_directories(${JSONCPP_INCLUDE_DIRS})

# Add directory to the current static shared library
target_include_directories(inputsJson PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Linking jsonCpp library
# target_link_libraries(inputsJson jsoncpp)