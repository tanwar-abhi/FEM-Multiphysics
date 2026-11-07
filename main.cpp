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
#include "material.hpp"

// #include <spdlog/spdlog.h>
// Incase we wish to export the logs to a file we can use the sink/basic_file_sink.h
// #include <spdlog/sinks/basic_file_sink.h>

#include <iomanip>
#include <string>
#include <filesystem>
#include <cstdlib>

void setSpdLogLevel(const std::string& envLogLevel);

int main()
{

    setSpdLogLevel("LOG_LEVEL");

    spdlog::info("Starting Multiphysics Simulator");
    const std::string inputsPath = std::filesystem::current_path().string() + "/sampleInputs/";
    spdlog::debug("(main) inputsPath = {}", inputsPath);

    //Object containing mesh parts defined by the user in the parts Input file
    PartsInput partsInput;
    partsInput.readPartInputs(inputsPath);
    // partsInput.printPartInputs();

    // Oject containing solver inputs defined by user in the inputs file
    SolverInput solverInput;
    if (!solverInput.readInputs(inputsPath))
    {
        spdlog::critical("(main) solver Inputs parsing failed");
        return -2;
    }

    Material materialDetails(inputsPath);

    if (materialDetails.readMaterialInputs())
    {
        spdlog::info("(main) material configs read successfully from file = {}", materialDetails.getMaterialInputFileName());
    }
    else {
        spdlog::critical("(main) Failed to read materials data from file = {}", materialDetails.getMaterialInputFileName());
        return -2;
    }

    // // LinearElasticMaterial solid(inputsPath);
    // std::vector<Solid> materialsSolid(inputsPath);
    // std::vector<MaterialThermal> materialsThermal(inputsPath);


    // // Read Material properties based on the equations defined in user input files
    // for (int i = 0; i < solverInput.getTotalEquations(); i++)
    // {
    //     if (solverInput.equations[i]->solverEquation == EquationType::PlaneStress || 
    //         solverInput.equations[i]->solverEquation == EquationType::PlaneStrain ||
    //         solverInput.equations[i]->solverEquation == EquationType::LinearElastic3D ) {
    //             // Solid solid;
    //             Solid solid;
    //             materialsSolid = solid.readMaterialInput();
    //         }
    // }

    return 0;
}



void setSpdLogLevel(const std::string& envVariable)
{
    char* envLogLevel = std::getenv(envVariable.c_str());

    if (envLogLevel != nullptr)
    {
        std::string logLevel(envLogLevel);
        std::transform(logLevel.begin(), logLevel.end(), logLevel.begin(), ::toupper);
        std::cout<<"env variable logLevel defined = "<<logLevel<<std::endl;

        if (logLevel == "TRACE") {
            spdlog::set_level(spdlog::level::trace);
        }
        else if (logLevel == "DEBUG") {
            spdlog::set_level(spdlog::level::debug);
        }
        else {
            spdlog::set_level(spdlog::level::info);
        }
        std::cout<<"Setting Log level = "<<logLevel<<std::endl;
    }
    else{
        std::cout<<"LogLevel not defined, setting default logLevel = INFO\n";
        spdlog::set_level(spdlog::level::info);
    }
}