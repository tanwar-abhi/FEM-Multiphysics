// Data Structures for solver input files

#pragma once

#include <jsoncpp/json/json.h>
#include <jsoncpp/json/value.h>

#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>

class Equation {

    public:
    Equation();
    ~Equation();

    std::string name;

    /*
    solverEq = 1 -> Linear Elastic (Plane Stress), 2 -> heat transfer, 3-> beam, 4->viscous incompressible,
    5->invicid incompressible, 6->viscous compressible, 7->inviscid compressible, 8-> Energy equation,
    9-> plate, 10-> Truss, 11->Frames (beams+Truss), 12-> Shell, 13-> Linear Elastic (Plane Strain),
    14 -> 3D Linear Elasticity, 15 -> Topology Optimization
    */
    int solverEq;

    // Mesh file path and name with extension {.msh}
    std::string meshFile;
    std::string meshField;

    // Name of material property
    std::string materialPropName;

    // Degree of freedom at each node for each problem type (as per solverEq)
    int DOF;

    /* Type of element is same as the elementType defined in gmsh document i.e. 
    1D Elements = 2NodeRod {gmsh=1}, 3NodeLine {gmsh=8} , 2NodeBeam
    2D Elements = 3NodeTria {gmsh=2}, 4NodeQuad {gmsh=3}, 6NodeTria {gmsh=9}, 9NodeQuad {gmsh=10}
    3D Elements = 8NodeHexahedron {gmsh=5}
    */
    int ElementType;

    // No. of Gauss Points {integration points} for each element
    int NGP;
    
    // Elemental Tag Id
    int elemTagId;

    // volume fraction
    double volumeFraction = 0;


    // Optimality Criteria Type
    int ocType = 0;
};

class SolverInp {

    public:
    SolverInp();
    ~SolverInp();

    // Copy constructor
    SolverInp(const SolverInp &obj);

    std::string coordinateSystem;
    int dimension;


    // Boolean variable to check if problem is transient or steady state
    bool isTransient;

    // Type of algorithm for solving Linear system of equations 
    // 1 -> classical ; 2 -> Hybrid {Quantum + classical}; 3 -> Only Quantum
    int algorithm = 1;

    // Total Time {SI units seconds} for transient problems
    double StartTime = 0;
    double EndTime = 0;
    double TotalTime = 0;

    // Size of each Time step to be defined for transient problems {SI units seconds}
    double dt = 0;

    // Stopping criteria {eps}
    double eps = 0;

    /* Time discretization method selection
    THETA = 0 -> Forward Euler
    THETA = 1.0/2.0 -> Crank-Nickolson
    THETA = 2.0/3.0 -> Galerkin
    THETA = 1 -> Backward Euler
    */
    float THETA = 1;

    // Type of mass matrix
    // 0-> Lumped Mass matrix ; 1-> Consistent Mass Matrix
    int massMatrixType;

    // number of equations
    int nEquations;

    Equation* equations = NULL;

    void readInputs(); 
};


