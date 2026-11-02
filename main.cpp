// Main file for the implmentation of simulation of heat transfer in an axisymmetric sold domain.


// #include "neumanBC.hpp"
// #include "postProcess.hpp"
// #include "dirichletBC.hpp"
// #include "Matrix_Assemble.hpp"


// // #include <Eigen/IterativeLinearSolvers>
// #include<Eigen/SparseLU>
// // #include<Eigen/SparseCholesky>

#include "partInputs.hpp"
#include "solverInputs.hpp"
#include "solid.hpp"
#include "thermal.hpp"

#include <iomanip>
#include <string>
#include <filesystem>

int main()
{
    const std::string inputsPath = std::filesystem::current_path().string() + "/sampleInputs/";
    std::cout << "inputsPath = "<< inputsPath<<std::endl;

    //Object containing mesh parts defined by the user in the parts Input file
    PartsInput partsInput;
    partsInput.readPartInputs(inputsPath);
    // partsInput.printPartInputs();

    // Oject containing solver inputs defined by user in the inputs file
    SolverInput solverInput;
    if (!solverInput.readInputs(inputsPath))
    {
        std::cout<<"[ERROR] (main) solver Inputs parsing failed\n";
        return -2;
    }


    // LinearElasticMaterial solid(inputsPath);
    std::vector<Solid> materialsSolid(inputsPath);
    std::vector<MaterialThermal> materialsThermal(inputsPath);
    

    // Read Material properties based on the equations defined in user input files
    for (int i = 0; i < solverInput.getTotalEquations(); i++)
    {
        if (solverInput.equations[i]->solverEquation == EquationType::PlaneStress || 
            solverInput.equations[i]->solverEquation == EquationType::PlaneStrain ||
            solverInput.equations[i]->solverEquation == EquationType::LinearElastic3D ) {
                // Solid solid;
                Solid solid;
                materialsSolid = solid.readMaterialInput();
            }
    }




    return 0;
}

