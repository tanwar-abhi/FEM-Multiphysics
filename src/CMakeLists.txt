
set(mpsSRCCode
	bForce_Traction.cpp 	bForce_Traction.hpp
	boundary_HT.cpp 		boundary_HT.hpp
	boundaryInputs.cpp		boundaryInputs.hpp
	dirichletBC.cpp			dirichletBC.hpp
	material.cpp			material.hpp
	Matrix_Assemle.cpp		Matrix_Assemle.hpp
	neumanBC.cpp			neumanBC.hpp
	structureBoundary.cpp	structureBoundary.hpp
	Time_Discretization.cpp	Time_Discretization.hpp
)

# Create a static library (*.so) of element code for linking with final executable
add_library(srcFunction STATIC ${mpsSRCCode})

# Add directory to the current static shared library
target_include_directories(srcFunction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})