// Data Structures for solver input files

#pragma once

#include <jsoncpp/json/json.h>
#include <jsoncpp/json/value.h>

#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <memory>


// helper enum for finite element shapes
enum class ElementShapeType 
{
    TwoNodeRod,
    ThreeNodeTri,
    FourNodeQuad,
    FourNodeTetra,
    EightNodeHexa,
    ThreeNodeRod,
    SixNodeTri,
    NineNodeQuad,
    TenNodeTetra,
    TwentySevenNodeHexa,
    Unknown
};

enum class EquationType{
    Beam,
    Plate,
    Truss,
    Frame,
    Shell,
    PlaneStress,
    PlaneStrain,
    LinearElastic3D,
    HeatTransfer,
    TopologyOptimization,
    Axisymmetric2D,
    Unknown
};

enum class SolverType
{
    Steady,
    Transient,
    Unknown
};


class Equation 
{
    private:
    // Name of config file that contains the user input prarameters to be read.
    const std::string fileName = "solverConfig.json";

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
    EquationType solverEquation = EquationType::Unknown;

    // Mesh file path and name with extension {.msh}
    std::string meshFile;
    std::string meshField;

    // Name of material property
    std::string materialPropertyName;

    // Degree of freedom at each node for each problem type (as per solverEq)
    int DOF;

    /* Type of element is same as the elementType defined in gmsh document i.e. 
    1D Elements = 2NodeRod {gmsh=1}, 3NodeLine {gmsh=8} , 2NodeBeam
    2D Elements = 3NodeTria {gmsh=2}, 4NodeQuad {gmsh=3}, 6NodeTria {gmsh=9}, 9NodeQuad {gmsh=10}
    3D Elements = 8NodeHexahedron {gmsh=5}
    */
    ElementShapeType ElementType = ElementShapeType::Unknown;

    // No. of Gauss Points {integration points} for each element
    int numberOfGaussPoints{0};
    
    // Elemental Tag Id
    int elementTagId;


    // volume fraction for Topology Optimization
    double volumeFraction = 0;

    // Optimality Criteria Type
    int ocType = 0;

    int penalization {0};
    double filterRadius {0.0};
};

class SolverInput 
{
    private:
    bool parseEquationsFromSolverJson(const std::string& solverJsonFile);
    bool parseSolverFromSolverJson(const std::string& solverJsonFile);

    EquationType searchEquationType(const std::string& equationString);

    const std::string fileName = "solverConfig.json";

    // number of equations
    int totalEquations{0};

    // Boolean variable to check if problem is transient or steady state
    bool isTransient {false};

    public:
    SolverInput();
    ~SolverInput();

    // Copy constructor
    // SolverInput(const SolverInput &obj);

    std::string coordinateSystem;
    int dimension;

    // std::string solverType;
    SolverType solverType = SolverType::Unknown;

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
    float THETA = 1.0;

    // Type of mass matrix
    // 0-> Lumped Mass matrix ; 1-> Consistent Mass Matrix
    int massMatrixType;

    std::vector<std::shared_ptr<Equation>> equations;

    bool readInputs(const std::string& inputFilePath);

    int getTotalEquations() const;

    bool getIsTransientSolver() const;
};


