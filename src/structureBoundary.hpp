// The functions defined within this are to intialize Body (gravity) force
//  and impose Traction forces on the Neuman/Natural boundary of domain.

#pragma once
// #include "../Libraries/Eigen/Dense"
// #include "BasisFunction.hpp"
// #include "readMesh.hpp"
// #include <set>
// #include <cmath>
// #include "ProblemParameters.hpp"
// #include "Boundary_HT.hpp"
#include "Matrix_Assemble.hpp"


// Initialize a vector of Body Forces based on value defined by user
// inputs -> BForce : BodyForce vector; BFV : Body Force values
// void initializeBForce(Eigen::VectorXd &BForce, Parameters_LE input);


/* Imposing neuman boundary traction force as per user definition only for the Natural/Essential boundary node number
INPUTS ::
    fGlobal = Global force vector (double)
    meshElement = Index of meshElement array containing element data 
    mesh = mesh object containg mesh details and neuman boundary edge nodes
    solverObj = Object of SolverInp class
    boundary = Object contianing data related to boundary conditions as per boundary.JSON
    material = Object of LinearElasticMaterial class containing material property details
    eqn = iteration of equation, for reference
*/
// void imposeTraction(Eigen::VectorXd &fGlobal, Element meshElement, readMesh mesh, SolverInp solverObj, BoundaryConditions boundary, LinearElasticMaterial material, int eqn);
void imposeTraction(Eigen::VectorXd &fGlobal, const Element &meshElement, const readMesh &mesh, const SolverInp &solverObj, const NeumannBC &tractionBoundary, const LinearElasticMaterial &material, int eqn);
double temp_detJ(const Eigen::MatrixXd &Xe);

