
#include "solverInputs.hpp"


Equation::Equation(){
    // Default constructor for equation class.
}

Equation::~Equation(){
    // Default destructor
}

SolverInput::SolverInput(){
    // Default Constructor of Solver inputs class.
}

SolverInput::~SolverInput(){
    // Default Destructor

    // std::cout<<"SolverInput destructor called\n";
    // Delete the pointer to release the heap block of memory
    // delete[] equations;
}


// Deep Copy constructor for SolverInputs, called when one object of class is instantiated using another object of class.
// It's always called when object of class pased by value in a function
SolverInput::SolverInput(const SolverInput &obj)
{
    coordinateSystem = obj.coordinateSystem;
    dimension = obj.dimension;
    isTransient = obj.isTransient;
    algorithm = obj.algorithm;
    StartTime = obj.StartTime; EndTime = obj.EndTime; TotalTime = obj.TotalTime;
    dt = obj.dt;
    eps = obj.eps; THETA = obj.THETA;
    massMatrixType = obj.massMatrixType;
    nEquations = obj.nEquations;
    equations = obj.equations;
}



// helper enum for finite element shapes
enum class FE_ShapeType {
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

// mapping from FE shape string to enum; faster and centralised
const std::unordered_map<std::string, FE_ShapeType> shapeMap = {
    {"2NodeRod", FE_ShapeType::TwoNodeRod},
    {"3NodeTri", FE_ShapeType::ThreeNodeTri},
    {"4NodeQuad", FE_ShapeType::FourNodeQuad},
    {"4NodeTetra", FE_ShapeType::FourNodeTetra},
    {"8NodeHexa", FE_ShapeType::EightNodeHexa},
    {"3NodeRod", FE_ShapeType::ThreeNodeRod},
    {"6NodeTri", FE_ShapeType::SixNodeTri},
    {"9NodeQuad", FE_ShapeType::NineNodeQuad},
    {"10NodeTetra", FE_ShapeType::TenNodeTetra},
    {"27NodeHexa", FE_ShapeType::TwentySevenNodeHexa}
};

// convert string to enum value using lookup table
static FE_ShapeType shapeFromString(const std::string &s)
{
    auto it = shapeMap.find(s);
    if (it != shapeMap.end())
        return it->second;
    return FE_ShapeType::Unknown;
}

int getElementType(const std::string &FE_Shape)
{
    FE_ShapeType shape = shapeFromString(FE_Shape);

    switch (shape) {
    case FE_ShapeType::TwoNodeRod:       return 1;
    case FE_ShapeType::ThreeNodeTri:     return 2;
    case FE_ShapeType::FourNodeQuad:     return 3;
    case FE_ShapeType::FourNodeTetra:    return 4;
    case FE_ShapeType::EightNodeHexa:    return 5;
    case FE_ShapeType::ThreeNodeRod:     return 8;
    case FE_ShapeType::SixNodeTri:       return 9;
    case FE_ShapeType::NineNodeQuad:     return 10;
    case FE_ShapeType::TenNodeTetra:     return 11;
    case FE_ShapeType::TwentySevenNodeHexa: return 12;
    default:
        std::cerr << "Error : Element selection unavailable, check the selected element" << std::endl;
        exit(-404);
    }

    // unreachable, but keeps compiler happy
    return 0;
}


enum class EquationType{
    PlaneStress,
    HeatTransfer,
    Beam,
    Plate,
    Truss,
    Frame,
    Shell,
    PlaneStrain,
    LinearElastic3D,
    Unknown
};


int getEquationType(std::string tmp)
{
    if (tmp == "PLANESTRESS" || tmp == "ELASTIC" || tmp == "LINEARELASTIC" || tmp == "LINEAR_ELASTIC")
        // LINEAR ELASTIC plane stress Equation
        return 1;
    else if (tmp == "HT" || tmp == "HEATTRANSFER" || tmp == "HEAT_TRANSFER")
        // Thermal Equation, HEAT TRANSFER problems
        return 2;
    else if (tmp == "BEAM")
        // std::cout << "Equation is BEAM" << std::endl;
        return 3;
    else if (tmp == "PLATE")
        // std::cout<<"Equation is PLATE" <<std::endl;
        return 9;
    else if (tmp == "TRUSS")
        // std::cout<<"Equation is TRUSS" <<std::endl;
        return 10;
    else if (tmp == "FRAME")
        // std::cout<<"Equation is FRAME" <<std::endl;
        return 11;
    else if (tmp == "SHELL")
        // std::cout<<"Equation is SHELL" <<std::endl;
        return 12;
    else if (tmp == "PLANESTRAIN" || tmp == "PE" || tmp == "LE_PLAINSTRAIN")
        // std::cout<<"Equation is LE PLAIN strain" <<std::endl;
        return 13;
    else if (tmp == "3DLINEARELASTIC" || tmp == "3DLE" || tmp == "LINEARELASTIC3D")
        // 3D linear Elasticity
        return 14;
    // else if (tmp == "OPTIMIZATION")
        // return 15;
    else{
        std::cerr<<"Error : Equation type not defined, "<<std::endl;
        exit(-417);
    }
    return 0;
}


void SolverInput::readInputs(const std::string& inputFilePath)
{
    std::ifstream solverJsonFile(inputFilePath + "solver.json");

    if (!solverJsonFile.is_open()){
        std::cout << "ERROR Solver Inputs :: could not open solver input file" << std::endl;
        exit(-404);
    }

    Json::Value solverRoot;
    Json::Reader reader;
    reader.parse(solverJsonFile, solverRoot);

    const Json::Value inputEquation =  solverRoot["equation"];
    nEquations = inputEquation.size();

    // Temporary string variable to read inputs
    std::string tmp;

    for (int index = 0; index < nEquations; ++index)
    {
        std::shared_ptr<Equation> eqn = std::make_shared<Equation>();
        eqn->name = inputEquation[index]["name"].asString();

        tmp = inputEquation[index]["type"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        eqn->solverEquation = getEquationType(tmp);

        tmp = inputEquation[index]["meshField"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        eqn->meshField = tmp;

        eqn->meshFile = inputEquation[index]["meshFile"].asString();

        tmp = inputEquation[index]["materialPropertyName"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        eqn->materialPropertyName = tmp;

        eqn->numberOfGaussPoints = inputEquation[index]["gaussPoints"].asInt();

        eqn->ElementType = getElementType(inputEquation[index]["feShape"].asString());
        if (equations[index]->solverEquation == 15)
        {
            // eqn.penalization = inputEquation[index]["penalization"].asInt();
            // eqn.filterRadius = inputEquation[index]["filterRadius"].asDouble();
            eqn->volumeFraction = inputEquation[index]["volumeFraction"].asDouble();
            // eqn.ocType = inputEquation[index]["ocType"].asInt();
        }
    }

    const Json::Value inputSolver = solverRoot["solver"];

    if (inputSolver.size() > 1){
        std::cerr << "Solver Input Error : The number of Solver keywords in solver input should be 1" << std::endl;
        exit(-400);
    }

    //   Solver Type ###### ("Solver" Keyword in solver file)
    // index = 0, beacaus there is always single solver property for any problem
    tmp =  inputSolver["type"].asString();
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

    if ( tmp == "TRANSIENT" ) 
    {
        // std::cout<<"The solver is transient" <<std::endl;
        isTransient = true;

        StartTime = inputSolver["startTime"].asDouble();
        EndTime = inputSolver["endTime"].asDouble();
        dt = inputSolver["timeStep"].asDouble();
        TotalTime = EndTime - StartTime;
        eps = inputSolver["stoppingCriteria"].asDouble();
        massMatrixType = inputSolver["massMatrixType"].asInt();
    }
    else if (tmp == "STEADY" || tmp == "STEADYSTATE" || tmp == "STEADY_STATE")
    {
        std::cout << "The solver is steady state" << std::endl;
        isTransient = false;
    }

    //  Solver dimension
    tmp = inputSolver["coordinateSystem"].asString();
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
    coordinateSystem = tmp;

    if (tmp == "1D")
        dimension = 1;
    else if (tmp == "2D" || tmp == "AXIS")
        dimension = 2;
    else if (tmp == "3D")
        dimension = 3;
    else{
        std::cerr<<"Dimension Error : Invalid Coordinate System selected"<<std::endl;
        exit(-403);
    }

    // Algorithm Selection
    // 1 -> classical ; 2 -> Hybrid {Quantum + classical}; 3 -> Only Quantum
    tmp = inputSolver["algorithm"].asString();
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

    // Ternary operator to make sure algorithm type is classical by default.
    algorithm = tmp == "QUANTUM" ? 3 : 1;

    // if (tmp == "CLASSICAL")
    //     algorithm = 1;
    // else if(tmp == "HYBRID")
    //     algorithm = 2;
    // else if(tmp == "QUANTUM"){
    //     algorithm = 3;
    // }

    for (int i = 0; i < nEquations; i++)
    {
        if (dimension == 2){
            // Linear elastic plane stress(1) and plane strain(13)
            if (equations[i]->solverEquation == 1 || equations[i]->solverEquation == 13){
                equations[i]->DOF = 2;
            }
            // Thermal equations (scalar field problems)
            else if (equations[i]->solverEquation == 2){
                equations[i]->DOF = 1;
            }
        }
        else if (dimension == 3)
        {
            // 3D linear elasticity equation
            if (equations[i]->solverEquation == 14){
                equations[i]->DOF = 3;
            }
            // Thermal problem (scalar field equation problem)
            else if (equations[i]->solverEquation == 2){
                equations[i]->DOF = 1;
            }
        }
        // 1D problems
        else{
            // ######### Define DOF for 1D problems here #############
        }
    }

    std::cout << "done reading the solver input file" << std::endl;
}
