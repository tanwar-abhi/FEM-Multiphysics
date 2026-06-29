
set(boundaryCode
    boundaryBase.cpp    boundaryBase.hpp
    dirichletBC.cpp     dirichletBC.hpp
    neumannBC.cpp       neumannBC.hpp
    initialBC.cpp       initialBC.hpp
)

# create a static libary (*.so) of the boundary coditions ode that can be used for linking with add_executable
add_library(boundaryCondtions STATIC ${boundaryCode})


# Add directory to the current static shared library
target_include_directories(boundaryCondtions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})