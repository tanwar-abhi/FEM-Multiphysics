
#include "BoundaryInputs.hpp"


DirichletBC::DirichletBC(){
    // Defalut constructor for DirichletBC class
}

DirichletBC::~DirichletBC(){
    // Defalut destructor for DirichletBC class
}

NeumannBC::NeumannBC(){
    // Default constructor for NeumannBC class
}

NeumannBC::~NeumannBC(){
    // Default Destructor for NeumannBC class
}

InitialBC::InitialBC(){
    // Default Constructor for InitialBC class
}

InitialBC::~InitialBC(){
    // Default destructor for InitialBC class
}

BoundaryConditions::BoundaryConditions(){
    // Default constructor for the boundaries super class
}

// Copy constructor, will be called whenever function is called in a function
// Deep copy (Thread safe) is performed to make sure pointers are copied correctly and meory released by destructor
BoundaryConditions::BoundaryConditions(const BoundaryConditions &obj)
{

    // Copy all elements of dirichlet boundary class objects
    if (obj.nDBC != 0){
        // Allocate proper new memory location for objects to copy
        dirichlet = new DirichletBC[obj.nDBC];

        for (int i = 0; i < obj.nDBC; i++){
            dirichlet[i] = obj.dirichlet[i];
        }
    }

    // Copy all elements of neumann boundary class objects
    if (obj.nNBC != 0){
        neumann = new NeumannBC[obj.nNBC];

        for (int i = 0; i < obj.nNBC; i++){
            neumann[i] = obj.neumann[i];
        }
    }

    // Copy all elements of Initial boundary class objects
    if (obj.nIBC != 0){
        initial = new InitialBC[obj.nIBC];

        for (int i = 0; i < obj.nIBC; i++)
        {
            initial[i] = obj.initial[i];
        }
    }
}


BoundaryConditions::~BoundaryConditions(){
    // Default destructor of the super class
    // std::cout<<" Boundary Inputs destructor called\n";

    // Delete the pointer to release the heap block of memory
    if (this->nDBC != 0)
        delete[] this->dirichlet;
    if (this->nNBC != 0)
        delete[] this->neumann;
    if (this->nIBC != 0)
        delete[] this->initial;
}


//Function reads the type (used in dirichlet and Neumann boundary conditions)
int getType(const Json::Value &inp_boundaryType, int index){

    //Variable for stroing boundary type entered by user
    std::string tmp;

    tmp = inp_boundaryType[index]["type"].asString();

    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

    if(tmp == "CONSTANT"){
        return 1;
    }
    else if(tmp == "TEMPORAL"){
        return 2;
    }
    else if(tmp == "CONSTANTNODALFILEINPUT"){
        return 3;
    }
    else if(tmp == "TEMPORALNODALFILEINPUT"){
        return 4;
    }
    else{
        std::cerr<<"INVALID TYPE :: Boundary type selected by user is invalid\n";
        exit(-403);
    }

    return 0;
}






//Function to read the Dirichlet boundary parameters from input file
void BoundaryConditions::readDirichletInputs()
{
    std::ifstream BC_text("boundary.json");
    Json::Reader BC_reader;
    Json::Value BC_root;

    bool parsingSuccessful = BC_reader.parse(BC_text, BC_root);

    if(!parsingSuccessful){
        std::cerr << "BOUNDARY INPUTS ERROR :: Error while parsing Dirichlet Boundary Conditions\n";
        exit(-404);
    }

    Json::Value inp_dirichletBC = BC_root["dirichletBC"];

    this->nDBC = inp_dirichletBC.size();
    DirichletBC *DBC = new DirichletBC[nDBC];

    for (int index = 0; index < nDBC; index++){

        std::string tmp = inp_dirichletBC[index]["name"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        DBC[index].name = tmp;

        DBC[index].type = getType(inp_dirichletBC, index);

        tmp = inp_dirichletBC[index]["variable"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        DBC[index].variable = tmp;

        // If body force defined in the boundary conditions JSON
        if (DBC[index].variable == "BODYFORCE"){
            for (int i = 0, n = inp_dirichletBC[index]["values"].size(); i < n; i++){
                DBC[index].values[i] = inp_dirichletBC[index]["values"][i].asDouble();
            }
            // bodyForceDefined = true;
        }
        else{
            DBC[index].value = inp_dirichletBC[index]["value"].asDouble();
        }

        DBC[index].meshFile = inp_dirichletBC[index]["meshFile"].asString();

        tmp = inp_dirichletBC[index]["meshField"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        DBC[index].meshField = tmp;
    }
    this->dirichlet = DBC;
}




// Function to read inital boundary conditions
void BoundaryConditions::readInitialInputs()
{
    std::ifstream BC_text("boundary.json");
    Json::Reader BC_reader;
    Json::Value BC_root;

    bool parsingSuccessful = BC_reader.parse(BC_text, BC_root);

    if(!parsingSuccessful){
        std::cout << "Error while parsing {IBC}\n";
        exit(-404);
    }

    Json::Value inp_initialBC = BC_root["initialBC"];
    std::string tmp;

    this->nIBC = inp_initialBC.size();

    InitialBC *IBC = new InitialBC[this->nIBC];

    for(int index = 0; index < nIBC; index++){

        IBC[index].name = inp_initialBC[index]["name"].asString();

        tmp = inp_initialBC[index]["type"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

        if (tmp == "CONSTANT")
            IBC[index].type = 1;
        else if (tmp == "NODALFILEINPUT"){
            IBC[index].type = 2;
        }

        tmp = inp_initialBC[index]["variable"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        IBC[index].variable = tmp;

        // If body force is defined in the problem
        if (IBC[index].variable == "BODYFORCE"){
            for (int i = 0, n = inp_initialBC[index]["values"].size(); i < n; i++){
                IBC[index].values[i] = inp_initialBC[index]["values"][i].asDouble();
            }
            // bodyForceDefined = true;
        }
        else{
            IBC[index].value = inp_initialBC[index]["value"].asDouble();
        }

        IBC[index].meshFile = inp_initialBC[index]["meshFile"].asString();
        tmp = inp_initialBC[index]["meshField"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        IBC[index].meshField = tmp;
    }

    this->initial = IBC;
}


// Read the Neumann boundary conditions defined in the json input file
void BoundaryConditions::readNeumannInputs()
{
    std::ifstream BC_text("boundary.json");
    Json::Reader BC_reader;
    Json::Value BC_root;

    //test to check if parsing is happening properly
    bool parsingSuccessful = BC_reader.parse(BC_text, BC_root);
    if(!parsingSuccessful){
        std::cout << "Error while parsing Neuman boundary inputs from json input file"<<std::endl;
        exit(-404);
    }

    Json::Value inp_neumannBC = BC_root["neumannBC"];

    this->nNBC = inp_neumannBC.size();

    //temporary variable for storing type of Problem (Solid Mechanics or Heat Transfer)
    std::string tmp;


    NeumannBC *NBarray = new NeumannBC[this->nNBC];

    //Loop for reading multiple inputs given by the user
    for(int index = 0; index < inp_neumannBC.size(); index++)
    {
        tmp = inp_neumannBC[index]["name"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        NBarray[index].name = tmp;

        NBarray[index].type = getType(inp_neumannBC, index);

        tmp = inp_neumannBC[index]["boundaryType"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        NBarray[index].boundaryType = tmp;

        tmp = inp_neumannBC[index]["variable"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        NBarray[index].variable = tmp;

        NBarray[index].meshFile = inp_neumannBC[index]["meshFile"].asString();
        tmp = inp_neumannBC[index]["meshField"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        NBarray[index].meshField = tmp;

        // if (NBarray[index].type == 1 && NBarray[index].boundaryType == "NORMALTOBOUNDARY"){
        if (NBarray[index].type == 1){

            if (NBarray[index].variable == "CONVECTIVEHEATTRANSFER"){
                NBarray[index].ambientTemp = inp_neumannBC[index]["ambientTemp"].asDouble();
                NBarray[index].H = inp_neumannBC[index]["convectionCoeff"].asDouble();
            }
            else if (NBarray[index].variable == "HEATSOURCESINK"){
                // Heat Source (+ve value) and Sink (-ve Value)
                NBarray[index].Q = inp_neumannBC[index]["Q"].asDouble();
                NBarray[index].domainSource = inp_neumannBC[index]["domainSource"].asBool();
            }
            // For dispX,Y,Z or velX,Y,Z or temperature 
            else{
                NBarray[index].value = inp_neumannBC[index]["value"].asDouble();
            }
        }


        if (NBarray[index].variable == "TRACTION" || NBarray[index].variable == "PRESSURE" || NBarray[index].variable == "HEATFLUX")
        {
            if (NBarray[index].boundaryType == "NORMALTOBOUNDARY")
            {
                NBarray[index].value = inp_neumannBC[index]["value"].asDouble();
            }
            else if(NBarray[index].boundaryType == "COMPONENTS")
            {
                if (inp_neumannBC[index]["values"].size() > 0){
                    for (int i = 0, n = inp_neumannBC[index]["values"].size(); i < n; i++){
                        NBarray[index].values[i] = inp_neumannBC[index]["values"][i].asDouble();
                    }
                }
                else{
                    NBarray[index].values[0] = inp_neumannBC[index]["value1"].asDouble();
                    NBarray[index].values[1] = inp_neumannBC[index]["value2"].asDouble();
                    NBarray[index].values[2] = inp_neumannBC[index]["value3"].asDouble();
                }
            }
            else
            {
                throw std::runtime_error("Heatflux/Pressure/Traction boundary type not defined");
            }
        }
        else if (NBarray[index].variable == "FORCE")
        {
            if(NBarray[index].boundaryType == "COMPONENTS")
            {
                if (inp_neumannBC[index]["values"].size() > 0){
                    for (int i = 0, n = inp_neumannBC[index]["values"].size(); i < n; i++){
                        NBarray[index].values[i] = inp_neumannBC[index]["values"][i].asDouble();
                    }
                }
                else{
                    NBarray[index].values[0] = inp_neumannBC[index]["value1"].asDouble();
                    NBarray[index].values[1] = inp_neumannBC[index]["value2"].asDouble();
                    NBarray[index].values[2] = inp_neumannBC[index]["value3"].asDouble();
                }
            }
            else
            {
                throw std::runtime_error("Force boundary type not defined");
            }
        }

    }

    this->neumann = NBarray;
}






// Function to read all boundary conditions defined in the boundary (json) input file
void BoundaryConditions::readBoundaryInputs()
{
    this->readDirichletInputs();
    this->readNeumannInputs();
    this->readInitialInputs();
}



