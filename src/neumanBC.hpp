
#pragma once
#include "boundary_HT.hpp"
// #include "Element.hpp"
// #include "SolverInputs.hpp"
#include "structureBoundary.hpp"


/* Apply all thermal neuman boundary conditions defined by user in json input files
Input Arguments : 
    tripletStiffness = Triplet (hashmap) containing identifiers (rows, columns) and values at those places.
    f_Global = R.H.S Vector containing thermal neuman boundary condition terms
    boundary = Object containing all boundary conditions and respective constants read from boundary.JSON
    meshElements[] = Array of meshElements, index of array represent the domain and boundary elements of domain
    mesh = Object of readMesh containing mesh data.
    solverInpObj = Object of SolverInp class containing solver related variables.
*/
void applyThermalNBC(std::vector<Tr> &tripletStiffness, Eigen::VectorXd &f_Global, const BoundaryConditions &boundary, 
                        const Element meshElements[], const readMesh &mesh, const SolverInp &solverInpObj, const MaterialThermal &material);




/* Apply all Structural neuman boundary conditions defined by user in problem
Input Arguments:
    f = R.H.S vector containing loads to impose
    mesh = Object of readMesh containing mesh data
    meshElement[] = Array of Elements of mesh, each index correspont to one of either domain or boundary elements
    solverInp = SolverInp class object containing data related to solver
    material_LE = Material properties of the domain
    boundary = Object containing all boundary conditions and respective constants read from boundary.JSON
    eqnIndex = Index number of the equation (in coupled problems)
*/
void applyStructureNBC(Eigen::VectorXd &f, const readMesh &mesh, const Element meshElements[], const SolverInp &solverInpObj, const LinearElasticMaterial &material_LE, const BoundaryConditions &boundary, int eqnIndex);


