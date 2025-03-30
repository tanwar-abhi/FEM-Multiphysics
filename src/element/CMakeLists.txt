
set(elementCode
	basisFunction.cpp	basisFunction.hpp
	element.cpp			element.hpp
	quadratures.cpp		quadratures.hpp
)

# Create a static library (*.so) of element code for linking with final executable
add_library(element STATIC ${elementCode})

# Add directory to the current static shared library
target_include_directories(element PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})