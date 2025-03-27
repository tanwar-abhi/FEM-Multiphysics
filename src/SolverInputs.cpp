
#include "SolverInputs.hpp"


Equation::Equation(){
    // Default constructor for equation class.
}

Equation::~Equation(){
    // Default destructor
}

SolverInp::SolverInp(){
    // Default Constructor of Solver inputs class.
}

SolverInp::~SolverInp(){
    // Default Destructor

    // std::cout<<"SolverInp destructor called\n";
    // Delete the pointer to release the heap block of memory
    delete[] equations;
}


// Deep Copy constructor for SolverInputs, called when one object of class is instantiated using another object of class.
// It's always called when object of class pased by value in a function
SolverInp::SolverInp(const SolverInp &obj)
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

    // Allocate new memory address to copy all equations
    equations = new Equation[nEquations];
    for (int i = 0; i < nEquations; i++)
    {
        equations[i] = obj.equations[i];
    }

}



int getElementType(std::string FE_Shape)
{
    if (FE_Shape == "2NodeRod")
        return 1;
    else if (FE_Shape == "3NodeTri")
        return 2;
    else if (FE_Shape == "4NodeQuad")
        return 3;
    else if (FE_Shape == "4NodeTetra")
        return 4;
    else if (FE_Shape == "8NodeHexa")
        return 5;
    else if (FE_Shape == "3NodeRod")
        return 8;
    else if (FE_Shape == "6NodeTri")
        return 9;
    else if (FE_Shape == "9NodeQuad")
        return 10;
    else if (FE_Shape == "10NodeTetra")
        return 11;
    else if (FE_Shape == "27NodeHexa")
        return 12;
    else{
        std::cerr<<"Error : Element selection unavailable, check the selected element"<<std::endl;
        exit(-404);
    }

    return 0;
}





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




void SolverInp::readInputs()
{

    // read  input files
    std::ifstream text("solver.json");

    if (!text.is_open()){
        std::cout << "ERROR Solver Inputs :: could not open solver input file" << std::endl;
        exit(-404);
    }

    Json::Value solver_root;
    Json::Reader reader;
    reader.parse(text, solver_root);

    Json::Value inp_equation =  solver_root["equation"];
    this->nEquations = inp_equation.size();

    // Temporary string variable to read inputs
    std::string tmp;

    Equation *eqn = new Equation [inp_equation.size()];
    this->equations = eqn;

    for (int index = 0; index < nEquations; ++index)
    {
        this->equations[index].name = inp_equation[index]["name"].asString();

        tmp = inp_equation[index]["type"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        this->equations[index].solverEq = getEquationType(tmp);

        tmp = inp_equation[index]["meshField"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        this->equations[index].meshField = tmp;

        this->equations[index].meshFile = inp_equation[index]["meshFile"].asString();

        tmp = inp_equation[index]["materialPropertyName"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        this->equations[index].materialPropName = tmp;

        this->equations[index].NGP = inp_equation[index]["gaussPoints"].asInt();

        this->equations[index].ElementType = getElementType(inp_equation[index]["feShape"].asString());
        if (equations[index].solverEq == 15)
        {
            // this->equations[index].penalization = inp_equation[index]["penalization"].asInt();
            // this->equations[index].filterRadius = inp_equation[index]["filterRadius"].asDouble();
            this->equations[index].volumeFraction = inp_equation[index]["volumeFraction"].asDouble();
            // this->equations[index].ocType = inp_equation[index]["ocType"].asInt();
        }
    }


    Json::Value inp_solver = solver_root["solver"];

    if (inp_solver.size() > 1){
        std::cerr << "Solver Input Error : The number of Solver keywords in solver input should be 1" << std::endl;
        exit(-400);
    }

    //   Solver Type ###### ("Solver" Keyword in solver file)
    // index = 0, beacaus there is always single solver property for any problem
    tmp =  inp_solver[0]["type"].asString();

    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

    if (tmp=="TRANSIENT") {
        // std::cout<<"The solver is transient" <<std::endl;
        this->isTransient = true;

        this->StartTime = inp_solver[0]["startTime"].asDouble();
        this->EndTime = inp_solver[0]["endTime"].asDouble();
        this->dt = inp_solver[0]["timeStep"].asDouble();
        this->TotalTime = this->EndTime - this->StartTime;
        this->eps = inp_solver[0]["stoppingCriteria"].asDouble();
        this->massMatrixType = inp_solver[0]["massMatrixType"].asInt();
    }
    else if (tmp == "STEADY" || tmp == "STEADYSTATE" || tmp == "STEADY_STATE")
    {
        // std::cout << "The solver is steady" << std::endl;
        this->isTransient = false;
    }

    // Algorithm Selection
    // 1 -> classical ; 2 -> Hybrid {Quantum + classical}; 3 -> Only Quantum
    tmp = inp_solver[0]["algorithm"].asString();
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

    //  Solver dimension
    tmp = inp_solver[0]["coordinateSystem"].asString();
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


    // Read algorithm type from the solver input file
    tmp = inp_solver[0]["algorithm"].asString();
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
            if (equations[i].solverEq == 1 || equations[i].solverEq == 13){
                equations[i].DOF = 2;
            }
            // Thermal equations (scalar field problems)
            else if (equations[i].solverEq == 2){
                equations[i].DOF = 1;
            }
        }
        else if (dimension == 3)
        {
            // 3D linear elasticity equation
            if (equations[i].solverEq == 14){
                equations[i].DOF = 3;
            }
            // Thermal problem (scalar field equation problem)
            else if (equations[i].solverEq == 2){
                equations[i].DOF = 1;
            }
        }
        // 1D problems
        else{
            // ######### Define DOF for 1D problems here #############
        }
    }

    // std::cout << "done reading the solver input file" << std::endl;
}

