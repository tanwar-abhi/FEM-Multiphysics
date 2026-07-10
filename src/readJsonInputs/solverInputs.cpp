
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


// // Deep Copy constructor for SolverInputs, called when one object of class is instantiated using another object of class.
// // It's always called when object of class pased by value in a function
// SolverInput::SolverInput(const SolverInput &obj)
// {
//     coordinateSystem = obj.coordinateSystem;
//     dimension = obj.dimension;
//     isTransient = obj.isTransient;
//     algorithm = obj.algorithm;
//     StartTime = obj.StartTime; EndTime = obj.EndTime; TotalTime = obj.TotalTime;
//     dt = obj.dt;
//     eps = obj.eps; THETA = obj.THETA;
//     massMatrixType = obj.massMatrixType;
//     nEquations = obj.nEquations;
//     equations = obj.equations;
// }

// mapping from FE shape string to enum; faster and centralised
const std::unordered_map<std::string, ElementShapeType> shapeMap = {
    {"2NODEROD", ElementShapeType::TwoNodeRod},
    {"3NODETRI", ElementShapeType::ThreeNodeTri},
    {"4NODEQUAD", ElementShapeType::FourNodeQuad},
    {"4NODETETRA", ElementShapeType::FourNodeTetra},
    {"8NODEHEXA", ElementShapeType::EightNodeHexa},
    {"3NODEROD", ElementShapeType::ThreeNodeRod},
    {"6NODETRI", ElementShapeType::SixNodeTri},
    {"9NODEQUAD", ElementShapeType::NineNodeQuad},
    {"10NODETETRA", ElementShapeType::TenNodeTetra},
    {"27NODEHEXA", ElementShapeType::TwentySevenNodeHexa},
    {"UNKNOWN", ElementShapeType::Unknown}
};

// const std::unordered_map<std::string, ElementShapeType> shapeMap = {
const std::unordered_map<std::string, EquationType> equationTypeMap = {
    {"PLANESTRESS", EquationType::PlaneStress},
    {"ELASTIC", EquationType::PlaneStress},
    {"LINEARELASTIC", EquationType::PlaneStress},
    {"HEATTRANSFER", EquationType::HeatTransfer},
    {"BEAM", EquationType::Beam},
    {"PLATE", EquationType::Plate},
    {"TRUSS", EquationType::Truss},
    {"FRAME", EquationType::Frame},
    {"SHELL", EquationType::Shell},
    {"PLANESTRAIN", EquationType::PlaneStrain},
    {"LINEARELASTIC3D", EquationType::LinearElastic3D},
    {"OPTIMIZATION", EquationType::TopologyOptimization},
    {"UNKNONW", EquationType::Unknown}
};

// convert shape function from string to enum value using lookup table
static ElementShapeType lookUpElementShapeType(const std::string &s)
{
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    auto it = shapeMap.find(s);
    if (it != shapeMap.end())
        return it->second;
    return ElementShapeType::Unknown;
}

// convert equation type from string to enum value using lookup table
EquationType searchEquationType(const std::string& equationString)
{
    std::transform(equationString.begin(), equationString.end(), equationString.begin(), ::toupper);
    auto it = equationTypeMap.find(equationString);
    if (it != equationTypeMap.end()) {
        return it->second;
    }
    return EquationType::Unknown;
}

bool SolverInput::readInputs(const std::string& inputFilePath)
{
    bool result{false};

    try{
        result = parseEquationsFromSolverJson(inputFilePath);
        if (result)
        {
            result = parseSolverFromSolverJson(inputFilePath);
            if (!result)
            {
                std::cout<<"[Error] Failed to parse solver details from sovlerInputs json file.\n";
                return result;
            }
            std::cout << "done reading the solver input file" << std::endl;
        }
        else
        {
            std::cout<<"[Error] Failed to parse Equations details from sovlerInputs json file.\n";
            return result;
        }
    }
    catch (const std::runtime_error& e)
    {
        result = false;
        std::cerr<<"[Error] Caught Runtime exception : "<<e.what()<<std::endl;
    }
    catch (const std::exception& e)
    {
        result = false;
        std::cerr<<"[Error] Caught Standard Exception : "<<e.what()<<std::endl;
    }
    catch (...) {
        result = false;
        std::cerr << "[Error] Caught an unknown exception!" << std::endl;
    }

    return result;
}


bool SolverInput::parseEquationsFromSolverJson(const std::string& jsonFileFullPath)
{
    std::ifstream solverJsonFile(jsonFileFullPath + fileName);

    std::cout<<"[Debug] SolverInput::parseEquationFromSolverJson path solverJsonFile = "<<jsonFileFullPath + fileName<<std::endl;

    if (!solverJsonFile.is_open()){
        std::runtime_error("[Error] SolverInput::parseEquationsFromSolverJson :: could not open solver input file");
    }

    Json::Value solverInput;
    Json::Reader reader;
    reader.parse(solverJsonFile, solverInput);

    const Json::Value equationInput = solverInput["equation"];
    totalEquations = equationInput.size();

    // Temporary string variable to read inputs
    std::string tmp;
    bool equationReadingResult {false};

    if (equationInput.isArray())
    {
        for (int index = 0; index < totalEquations; index++)
        {
            std::shared_ptr<Equation> eqn = std::make_shared<Equation>();

            if (equationInput[index].isMember("name") && equationInput[index]["name"].isString())
            {
                eqn->name = equationInput[index]["name"].asString();
            }
    
            if (equationInput[index].isMember("type") && equationInput[index]["type"].isString())
            {
                tmp = equationInput[index]["type"].asString();
                eqn->solverEquation = searchEquationType(tmp);
            }
            else{
                throw std::runtime_error("sovler input file, equation object doesn't contain equation Type.");
            }


            if (equationInput[index].isMember("meshField") && equationInput[index]["meshField"].isString())
            {
                tmp = equationInput[index]["meshField"].asString();
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
                eqn->meshField = tmp;
            }
            else{
                throw std::runtime_error("sovler input file, equation field doesn't contain meshField.");
            }
    
            if (equationInput[index].isMember("meshFile") && equationInput[index]["meshFile"].isString())
            {
                eqn->meshFile = equationInput[index]["meshFile"].asString();
            }
            else{
                throw std::runtime_error("sovler input file, equation field doesn't contain mesh file name.");
            }
    
            if (equationInput[index].isMember("materialPropertyName") && equationInput[index]["materialPropertyName"].isString())
            {
                tmp = equationInput[index]["materialPropertyName"].asString();
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
                eqn->materialPropertyName = tmp;
            }
            else{
                throw std::runtime_error("sovler input file, equation field doesn't contain materialPropertyName.");
            }

            if (equationInput[index].isMember("gaussPoints") && equationInput[index]["gaussPoints"].isInt())
            {
                eqn->numberOfGaussPoints = equationInput[index]["gaussPoints"].asInt();
            }
            else{
                throw std::runtime_error("sovler input file, equation field doesn't contain total gauss points for each element.");
            }

            eqn->ElementType = lookUpElementShapeType(equationInput[index]["feShape"].asString());

            if (eqn->solverEquation == EquationType::TopologyOptimization)
            {
                if (equationInput[index].isMember("volumeFraction") && equationInput[index]["volumeFraction"].isString())
                {
                    eqn->volumeFraction = equationInput[index]["volumeFraction"].asDouble();
                }
                if (equationInput[index].isMember("penalization") && equationInput[index]["penalization"].isString())
                {
                    eqn->penalization = equationInput[index]["penalization"].asInt();
                }
                if (equationInput[index].isMember("filterRadius") && equationInput[index]["filterRadius"].isString())
                {
                    eqn->filterRadius = equationInput[index]["filterRadius"].asDouble();
                }
                if (equationInput[index].isMember("ocType") && equationInput[index]["ocType"].isString())
                {
                    eqn->ocType = equationInput[index]["ocType"].asInt();
                }
            }
            std::cout<<equations.size()<<"\n";
            equations.emplace_back(eqn);
            // std::cout<<"after emplace_back"<<equations.size()<<"\n";x
            equationReadingResult = true;
        }
    }
    else{
        throw std::runtime_error("solver.json input file does't contain equation details.");
    }

    return equationReadingResult;
}


bool SolverInput::parseSolverFromSolverJson(const std::string& jsonFileFullPath)
{
    std::ifstream solverJsonFile(jsonFileFullPath + fileName);

    std::cout<<"[Debug] SolverInput::parseSolverFromSolverJson path solverJsonFile = "<<jsonFileFullPath + fileName<<std::endl;

    if (!solverJsonFile.is_open()){
        std::runtime_error("[Error] SolverInput::parseSolverFromSolverJson :: could not open solver input file");
    }
    
    Json::Value solverInput;
    Json::Reader reader;
    reader.parse(solverJsonFile, solverInput);

    const Json::Value inputSolver = solverInput["solver"];

    bool solverParsingSucessful {false};

    std::string tmp;

    if (!inputSolver.isObject()){
        std::runtime_error("[Error] SolverInput::parseSolverFromSovlerJson solver.json file does't contian required solver object.");
    }
    // if (inputSolver.isObject())
    else
    {
        //   Solver Type ###### ("Solver" Keyword in solver file)
        // index = 0, beacaus there is always single solver property for any problem
        if (inputSolver.isMember("type") && inputSolver["type"].isString())
        {
            tmp =  inputSolver["type"].asString();
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            tmp == "TRANSIENT" ? solverType = SolverType::Transient : solverType = SolverType::Steady;
        }
        else{
            throw std::runtime_error("Solver type is not present");
        }

        if (solverType == SolverType::Transient)
        {
            std::cout<<"The solver type is transient" <<std::endl;
            isTransient = true;
    
            if (inputSolver.isMember("startTime") && inputSolver["startTime"].isDouble())
            {
                StartTime = inputSolver["startTime"].asDouble();
            }
            if (inputSolver.isMember("endTime") && inputSolver.isDouble())
            {
                EndTime = inputSolver["endTime"].asDouble();
            }
            if (inputSolver.isMember("timeStep") && inputSolver.isDouble())
            {
                dt = inputSolver["timeStep"].asDouble();
            }
            if(inputSolver.isMember("massMatrixType") && inputSolver["massMatrixType"].isInt())
            {
                massMatrixType = inputSolver["massMatrixType"].asInt();
            }
            if(inputSolver.isMember("stoppingCriteria") && inputSolver["stoppingCriteria"].isDouble())
            {
                eps = inputSolver["stoppingCriteria"].asDouble();
            }
            TotalTime = EndTime - StartTime;
        }
        // else if (tmp == "STEADY" || tmp == "STEADYSTATE" || tmp == "STEADY_STATE")
        else if (solverType == SolverType::Steady)
        {
            std::cout << "The solver is steady state" << std::endl;
            isTransient = false;
        }
        else {
            throw std::runtime_error("[Error] Incorrect Solver type in solverConfig.json, currently only steady and transient solver available.");
        }
    
        //  Solver dimension
        if (inputSolver.isMember("coordinateSystem") && inputSolver["coordinateSystem"].isString())
        {
            tmp = inputSolver["coordinateSystem"].asString();
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            coordinateSystem = tmp;
        }
    }

    if (coordinateSystem == "1D")
        dimension = 1;
    else if (coordinateSystem == "2D" || coordinateSystem == "AXIS")
        dimension = 2;
    else if (coordinateSystem == "3D")
        dimension = 3;
    else{
        std::runtime_error("[Error] SovlerInput::parseSolverFromSolverJson Dimension Invalid Coordinate System selected");
    }

    std::cout<<"CoordinateSystem ="<<coordinateSystem<<" ; Dimension ="<<dimension<<std::endl;
    // // Algorithm Selection
    // // 1 -> classical ; 2 -> Hybrid {Quantum + classical}; 3 -> Only Quantum
    // tmp = inputSolver["algorithm"].asString();
    // std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

    // // Ternary operator to make sure algorithm type is classical by default.
    // algorithm = tmp == "QUANTUM" ? 3 : 1;

    // // if (tmp == "CLASSICAL")
    // //     algorithm = 1;
    // // else if(tmp == "HYBRID")
    // //     algorithm = 2;
    // // else if(tmp == "QUANTUM"){
    // //     algorithm = 3;
    // // }

    solverParsingSucessful = true;

    for (int i = 0; i < totalEquations; i++)
    {
        if (dimension == 2){
            // Linear elastic plane stress(1) and plane strain(13)
            if (equations[i]->solverEquation == EquationType::PlaneStress || equations[i]->solverEquation == EquationType::PlaneStrain){
                equations[i]->DOF = 2;
            }
            // Thermal equations (scalar field problems)
            else if (equations[i]->solverEquation == EquationType::HeatTransfer){
                equations[i]->DOF = 1;
            }
        }
        else if (dimension == 3)
        {
            // 3D linear elasticity equation
            if (equations[i]->solverEquation == EquationType::LinearElastic3D){
                equations[i]->DOF = 3;
            }
            // Thermal problem (scalar field equation problem)
            else if (equations[i]->solverEquation == EquationType::HeatTransfer){
                equations[i]->DOF = 1;
            }
        }
        // 1D problems
        else{
            // ######### Define DOF for 1D problems here #############
            solverParsingSucessful = false;
            throw std::runtime_error("1D problems are currently unavailable.");
        }
    }
    return solverParsingSucessful;
}

int SolverInput::getTotalEquations() const
{
    return totalEquations;
}

bool SolverInput::getIsTransientSolver() const
{
    return isTransient;
}