// Main file for the implmentation of simulation of heat transfer in an axisymmetric sold domain.


// #include "src/Time_Discretization.hpp"
// #include <fstream>
// #include "src/HeatTransfer.hpp"
// #include "src/SolverInputs.hpp"
// #include "src/Element.hpp"
// #include "src/BoundaryInputs.hpp"
// #include "src/PostProcess.hpp"
// #include "src/neumanBC.hpp"
// #include "src/PostProcess.hpp"
// #include "src/DirichletBC.hpp"
// #include "src/Matrix_Assemble.hpp"

#include "PostProcess.hpp"
#include "neumanBC.hpp"
#include "PostProcess.hpp"
#include "DirichletBC.hpp"
#include "Matrix_Assemble.hpp"


// #include <Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
// #include<Eigen/SparseCholesky>

#include <iomanip>
#include <string>

int main()
{

    //Object containing mesh parts defined by the user in the parts Input file
    PartsInp partsInpObj;
    partsInpObj.readPartInputs();

    // Oject containing solver inputs defined by user in the inputs file
    SolverInp solverInpObj;
    solverInpObj.readInputs();

    // Material vector objects containing material details for problem
    std::vector<MaterialThermal> materialThermal;
    std::vector<LinearElasticMaterial> material_LE;

    // Read Material properties based on the equations defined in user input files
    for (int i = 0; i < solverInpObj.nEquations; i++)
    {
        // Linearly elastic material
        if (solverInpObj.equations[i].solverEq == 1 || solverInpObj.equations[i].solverEq == 13 || solverInpObj.equations[i].solverEq == 14)
            material_LE = LinearElasticMaterial::readMaterialInputs();
        // Thermal problem
        else if (solverInpObj.equations[i].solverEq == 2){
            materialThermal = MaterialThermal::readMaterialInputs();
        }
    }


    // Create an object of super class containing all boundary conditions
    BoundaryConditions boundary;
    boundary.readBoundaryInputs();


    // Read mesh object constructor.
    readMesh mesh(partsInpObj, solverInpObj.dimension);

    // Read elements
    Element meshElements[mesh.NElementalTags];

    // Reading and Mapping mesh Elements with element Tag id's
    for (int tagCounter = 0; tagCounter < mesh.NElementalTags; tagCounter++)
    {
        meshElements[tagCounter].readElement(mesh, tagCounter);
    }

    // Mapping mesh with the equation
    for (int eqn = 0; eqn < solverInpObj.nEquations; eqn++)
    {
        for (int elemTagId = 0; elemTagId < mesh.NElementalTags; elemTagId++)
        {
            // std::cout<< solverInpObj.equations[eqn].meshFile << " " << mesh.elementalTagsMeshName[elemTagId] << " " << solverInpObj.equations[eqn].meshField << " " << mesh.elementalTags[elemTagId] << std::endl;
            if ((solverInpObj.equations[eqn].meshFile == mesh.elementalTagsMeshName[elemTagId]) && (solverInpObj.equations[eqn].meshField == mesh.elementalTags[elemTagId]))
            {
                solverInpObj.equations[eqn].elemTagId = elemTagId;
                break;
            }
        }
    }

    // Mapping mesh with the dirichlet boundaries
    for (int index = 0; index < boundary.nDBC; index++)
    {
        for (int elemTagId = 0; elemTagId < mesh.NElementalTags; elemTagId++)
        {
            if ((boundary.dirichlet[index].meshFile == mesh.elementalTagsMeshName[elemTagId]) && (boundary.dirichlet[index].meshField == mesh.elementalTags[elemTagId]))
            {
                //   std::cout<<"Mapping "<<mesh.elementalTags[elemTagId]<< " to the tagId "<<elemTagId<<std::endl;
                boundary.dirichlet[index].elemTagId = elemTagId;
                break;
            }
        }
    }

    // Mapping mesh with the initial boundaries
    if (solverInpObj.isTransient){
        for (int index = 0; index < boundary.nIBC; index++)
        {
            for (int elemTagId = 0; elemTagId < mesh.NElementalTags; elemTagId++)
            {
                if ((boundary.initial[index].meshFile == mesh.elementalTagsMeshName[elemTagId]) && (boundary.initial[index].meshField == mesh.elementalTags[elemTagId]))
                {
                    boundary.initial[index].elemTagId = elemTagId;
                    break;
                }
            }
        }
    }


    // Mapping mesh with the neumann boundaries
    for (int index = 0; index < boundary.nNBC; index++)
    {
        for (int elemTagId = 0; elemTagId < mesh.NElementalTags; elemTagId++)
        {
            if ((boundary.neumann[index].meshFile == mesh.elementalTagsMeshName[elemTagId]) && (boundary.neumann[index].meshField == mesh.elementalTags[elemTagId]))
            {
                boundary.neumann[index].elemTagId = elemTagId;
                break;
            }
        }
    }


    // Mapping material property of domain to equation
    for (int index = 0; index < solverInpObj.nEquations; index++){
        int solverEquation = solverInpObj.equations[index].solverEq;
        if (solverEquation == 2){
            for (int matIndex = 0, nMaterials = materialThermal.size(); matIndex < nMaterials; matIndex++){
                if (solverInpObj.equations[index].materialPropName == materialThermal[matIndex].name)
                    materialThermal[matIndex].elemTagId = solverInpObj.equations[index].elemTagId;
            }
        }
        else if (solverEquation == 1 || solverEquation == 13 || solverEquation == 14){
            for (int matIndex = 0, nMaterials = material_LE.size(); matIndex < nMaterials; matIndex++){
                if (solverInpObj.equations[index].materialPropName == material_LE[matIndex].name)
                    material_LE[matIndex].elemTagId = solverInpObj.equations[index].elemTagId;
            }
        }
    }


    // Value of Unknown at each node over the domain
    Eigen::MatrixXd Solution;

    // Global Stiffness Matrix and Force vector
    Eigen::SparseMatrix<double> K(mesh.NNodes, mesh.NNodes);
    // Eigen::SparseMatrix<double> K(mesh.NNodes*solverInpObj.equations[0].DOF, mesh.NNodes*solverInpObj.equations[0].DOF);
    K.setZero();

    Eigen::VectorXd f(mesh.NNodes * solverInpObj.equations[0].DOF);
    f.setZero();

    // Transient Solver
    if (solverInpObj.isTransient)
    {
        /*
        ** TODO be updated as per the new and imporved function
        */
    }

    // Steady State equation defined
    else
    {
        // dummy mass matrix (not required in steady state case, but defined as dummy variable to maintain functionality)
        // SpMatDouble C;

        // Vector of triplets containing only non zero element values and row, column indexes for global stiffness matrix
        std::vector<Tr> tripletStiffness;
        std::vector<Tr> tripletMass;


        for (int eqn=0; eqn < solverInpObj.nEquations; eqn++)
        {
            // Tag for domain element
            int domainElemTag = solverInpObj.equations[eqn].elemTagId;

            // Thermal problem stiffness matrix
            if (solverInpObj.equations[eqn].solverEq == 2)
            {
                // f = Eigen::VectorXd::Zero(mesh.NNodes * solverInpObj.equations[eqn].DOF);

                // Assemble the global stiffness matrix and Force vector for thermal problem using triplets
                assembleMatrixHT(tripletStiffness, tripletMass, f, mesh.Node_Coord, meshElements, domainElemTag, solverInpObj, materialThermal, boundary);

                // Apply thermal neumann Boundary conditions
                applyThermalNBC(tripletStiffness, f, boundary, meshElements, mesh, solverInpObj,  materialThermal[0]);

            }

            // Linear elastic problem striffness matrix
            if (solverInpObj.equations[eqn].solverEq == 1 || solverInpObj.equations[eqn].solverEq == 13 || solverInpObj.equations[eqn].solverEq == 14)
            {
                // f = Eigen::VectorXd::Zero(mesh.NNodes * solverInpObj.equations[eqn].DOF);

                assembleMatrixLE(tripletStiffness, f, mesh.Node_Coord, meshElements[domainElemTag], domainElemTag, solverInpObj, material_LE, boundary, eqn);

                // Apply structure neumann Boundary conditions
                applyStructureNBC(f, mesh, meshElements, solverInpObj, material_LE[0], boundary, eqn);

            }
        }

        std::cout<<"\nSteady state K matrix assembly Compleated "<< std::endl;


        // Converting triplet to sparse matrix (adding to sparse matrix)
        // K.setFromTriplets(tripletStiffness.begin(), tripletStiffness.end());


        // Apply all dirichlet boundary conditions defined by user in boundary.json and generate reduced stiffness matrix
        K = applyDirichletBC(tripletStiffness, f, meshElements, boundary, solverInpObj, mesh.NNodes);


        // ########### Check the time taken for compressed vs uncompressd CCS matrix ###########
        K.makeCompressed();

        // Solution of reduced system of linear equations
        Eigen::VectorXd Sol_Reduced;

        // ExportCSV(K, "results/A-thermal.csv");
        // ExportCSV(f, "results/f-thermal.csv");

        // Classical algorithm
        if (solverInpObj.algorithm == 1)
        {
            // Solver Selection
            // CG was used due to symmetric and positive definate stiffness matrix.
            // Eigen::ConjugateGradient<SpMatDouble> solver;
            Eigen::SparseLU<SpMatDouble> solver;
            // Eigen::BiCGSTAB<SpMatDouble, Eigen::IncompleteLUT<SpMatDouble>> solver;
            // Eigen::SimplicialCholesky<SpMatDouble> solver;
            // Eigen::SimplicialCholeskyLDLT<SpMatDouble> solver;
            // Eigen::SimplicialCholeskyLLT<SpMatDouble> solver;


            solver.compute(K);

            if(solver.info() != Eigen::ComputationInfo::Success){
                std::cout<<"solving failed"<<std::endl;
                return -408;
            }

            // Solve linear system of equation using ConjugateGradient solver
            Sol_Reduced = solver.solve(f);

            // ExportCSV(Sol_Reduced, "results/x-vector.csv");
        }
        // Qunatum algorithm
        else if (solverInpObj.algorithm == 3)
        {
            // Terminate the code to call upon the quantum linear solver
            // system("python quantumLinearSolver.py");
            return 101;
        }


        // Final solution after imposing value at Dirichlet boundary nodes.
        Eigen::MatrixXd Solution;

        if (boundary.nDBC != 0)
        {
            // Impose the dirichlet values to get the final solution
            Solution = imposeDirichletValues(Sol_Reduced, meshElements, boundary, mesh.NNodes, solverInpObj.equations[0].DOF);
        }
        else{
            Solution = Sol_Reduced;
        }

        // Output the solution to a file
        // ExportCSV(Solution, "results/Results.csv");


        // Post Processing for visualisation
        for (int eqn = 0; eqn < solverInpObj.nEquations; eqn++)
        {
            int domainElemTag = solverInpObj.equations[eqn].elemTagId;
            std::string resultsFile;

            // Thermal problem results
            if (solverInpObj.equations[eqn].solverEq == 2)
                resultsFile = "results/thermal.vtk";

            // Linear elastic problem results
            else if (solverInpObj.equations[eqn].solverEq == 1 || solverInpObj.equations[eqn].solverEq == 13 || solverInpObj.equations[eqn].solverEq == 14){
                resultsFile = "results/displacement.vtk";
            }
            GenerateVTK(resultsFile, mesh.Node_Coord, meshElements[domainElemTag], solverInpObj.equations[eqn], Solution);
        }
    }

    return 0;
}

