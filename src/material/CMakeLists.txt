
set(materialsDefinition
    MaterialBase.cpp    MaterialBase.hpp
    solid.cpp           solid.hpp
    thermal.cpp         thermal.hpp
)


# Create a static library (*.so) of element code for linking with final executable
add_library(materialsProperty STATIC ${materialsDefinition})

# Add directory to the current static shared library
target_include_directories(materialsProperty PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(materialsProperty PUBLIC ${CMAKE_SOURCE_DIR}/Libraries/eigen3.3)

