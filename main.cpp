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

#include <iomanip>
#include <string>
#include <filesystem>

int main()
{
    const std::string inputsPath = std::filesystem::current_path().string() + "/inputs/";

    //Object containing mesh parts defined by the user in the parts Input file
    PartsInput partsInput;
    partsInput.readPartInputs(inputsPath);

    // Oject containing solver inputs defined by user in the inputs file
    SolverInput solverInput;
    solverInput.readInputs(inputsPath);

    std::cout<<solverInput.dimension<<" "<<solverInput.equations[0]->meshField<<std::endl;


    return 0;
}

