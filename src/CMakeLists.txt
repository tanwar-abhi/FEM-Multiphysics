
set(mpsSRCCode
	boundary_HT.cpp 		boundary_HT.hpp
	boundaryInputs.cpp		boundaryInputs.hpp
	dirichletBC.cpp			dirichletBC.hpp
	material.cpp			material.hpp
	Matrix_Assemble.cpp		Matrix_Assemble.hpp
	neumanBC.cpp			neumanBC.hpp
	structureBoundary.cpp	structureBoundary.hpp
)

# Create a static library (*.so) of element code for linking with final executable
add_library(srcFunction STATIC ${mpsSRCCode})

add_subdirectory(element)
add_subdirectory(postProcessing)
add_subdirectory(preProcessing)
add_subdirectory(readJsonInputs)
 
# Add directory to the current static shared library
target_include_directories(srcFunction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(srcFunction PUBLIC ${CMAKE_SOURCE_DIR}/Libraries/eigen3.3)

# Linking jsonCpp library
target_link_libraries(srcFunction jsoncpp)
target_link_libraries(srcFunction element)