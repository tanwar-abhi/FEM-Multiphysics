
#pragma once

#include <jsoncpp/json/json.h>
#include <jsoncpp/json/value.h>

#include <vector>
// #include <string>
// #include "ProblemParameters.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>


enum class BoundaryType
{
    Constant,
    Temporal,
    ConstantNodalFileInput,
    TemporalNodalFileInput,
    Unknown
};



class BoundaryBase
{
    private:
    // Type of boundary Condition
    // 1=constant, 2=temporal, 3 = constant Nodal FileInput, 4 = temporal Nodal FileInput
    // int type {0};


    // variable defined on this [dirichlet] boundary
    // This parameter can take “dispX”, “dispY”, “dispZ”, “velX”, “velY”, “velZ”, "temperature", "bodyForce", "fixed", "free"
    std::string variable;

    // // Variable (or constants) defined on that particular boundary [Neumann]
    // // "traction", "convectiveHeatTransfer", "heatFlux", "HeatSourceSink", "dispX", "dispY", "dispZ", "velX", "velY", "velZ", "temperature"
    // std::string variable;

    // Elemental Tag Id
    int elementTagId;

    // name is a string, that will help to distinguish between multiple boundary conditions
    std::string name;

    public:
    BoundaryBase();
    ~BoundaryBase();

    // Path and Name of the mesh file
    std::string meshFile;

    // This is the physicalId defined in gmsh
    std::string meshField;

    BoundaryType getBoundaryType();
};