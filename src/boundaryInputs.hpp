
#pragma once

#include <jsoncpp/json/json.h>
#include <jsoncpp/json/value.h>

#include <vector>
// #include <string>
// #include "ProblemParameters.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>



class DirichletBC
{

    public:
    // Default Constructor and Destructor
    DirichletBC();
    ~DirichletBC();

    // Type of Dirichlet boundary
    // 1=constant, 2=temporal, 3 = constant Nodal FileInput, 4 = temporal Nodal FileInput
    int type;

    // variable defined on this dirichlet boundary
    // This parameter can take “dispX”, “dispY”, “dispZ”, “velX”, “velY”, “velZ”, "temperature", "bodyForce", "fixed", "free"
    std::string variable;

    // name is a string, that will help to distinguish between multiple boundary conditions
    std::string name;

    // Path and Name of the mesh file
    std::string meshFile;

    // This is the physicalId defined in gmsh
    std::string meshField;

    // Scalor value (depend upon type of variables)
    double value = 0;

    // array values (depend upon type of variables)
    double values[3] = {0,0,0};

    // Elemental Tag Id
    int elemTagId;

    // filename ????
    std::string filename;
};


class NeumannBC
{
    public:
    // Default Constructor and Destructor
    NeumannBC();
    ~NeumannBC();

    // Type of boundary condition
    // 1=constant, 2=temporal, 3=constant Nodal FileInput, 4=temporal Nodal FileInput
    int type;

    // Variable (or constants) defined on that particular boundary
    // "traction", "convectiveHeatTransfer", "heatFlux", "HeatSourceSink", "dispX", "dispY", "dispZ", "velX", "velY", "velZ", "temperature"
    std::string variable;

    //  This parameter can take "normalToBoundary", "gradX", "gradY", "gradZ"
    std::string boundaryType;

   // Name variable to identify the type of boundary condition
    std::string name;

    // Path and Name of the mesh file
    std::string meshFile;

    // This is the physicalId defined in gmsh
    std::string meshField;

    // Scalor value (depend upon type of variable)
    double value = 0;

    // Array of values (depend upon type of variable)
    double values[3] = {0,0,0};

    // Elemental Tag Id
    int elemTagId;

    // Convective heat transfer coefficient (H) {SI units in W/m2°C}
    double H = 0;

    // Convective heat exchange temperature i.e. ambient Temperature (ambientTemp) {SI units in °C}
    double ambientTemp = 0;

    // Internal heat generation rate per unit volume {SI units = W/m3}
    // positive(+) value of Q indicates an internal heat source; negative(-) value of Q means a heat sink (heat drawn out of volume)
    double Q = 0;

    // Flag to signify if heat source is over whole domain (true) or point/line source (false)
    bool domainSource;

    // filename ????
    std::string filename;
};



class InitialBC{
    public:
    InitialBC();
    ~InitialBC();

    // Boundary name is the physicalId defined in gmsh
    std::string name;

    // Initial value type
    // 1->constant, 2->Nodal File Input
    int type;

    // A string that contains the type of variable i.e. displacement, velocity , temperature
    std::string variable;

    // Path and Name of the mesh file
    std::string meshFile;

    // This is the physicalId defined in gmsh
    std::string meshField;

    // value (depend upon type of variables)
    double value = 0;

    // array of values (depend upon type of variables)
    double values[3] = {0,0,0};

    // Elemental Tag Id
    int elemTagId;

    // filename ????
    std::string filename;
};



class BoundaryConditions
{
    private:
    void readDirichletInputs();
    void readNeumannInputs();
    void readInitialInputs();

    public:
    BoundaryConditions();
    ~BoundaryConditions();

    // Copy constructor, will be called whenever function is called in a function
    BoundaryConditions(const BoundaryConditions &bcObj);


    DirichletBC *dirichlet = NULL;
    NeumannBC *neumann = NULL;
    InitialBC *initial = NULL;

    // Total number of total number of boundary conditions defined by the user.
    // nDBC-> Number of Dirichlet Boundary Condition; 
    int nDBC = 0;

    // nNBC-> Number of Neuman Boundary Condition; 
    int nNBC = 0;

    // nIBC-> Number of Initial Boundary Condition
    int nIBC = 0;

    void readBoundaryInputs();

    // Flag to identify if body force is defined in boundary conditions or not
    // bool bodyForceDefined = false;
};



