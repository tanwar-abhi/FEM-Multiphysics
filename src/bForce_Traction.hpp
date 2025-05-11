// The functions defined within this are to intialize Body (gravity) force
//  and impose Traction forces on the Neuman/Natural boundary of domain.

#pragma once

#include "basisFunction.hpp"
#include "readMesh.hpp"
#include "boundary_HT.hpp"
// #include "ProblemParameters.hpp"

// #include "../Libraries/Eigen/Dense"
#include <Eigen/Sparse>

#include <set>
#include <cmath>

// Initialize a vector of Body Forces based on value defined by user
// inputs -> BForce : BodyForce vector; BFV : Body Force values
void initializeBForce(Eigen::VectorXd &BForce, Parameters_LE input);


/* Imposing neuman boundary traction force as per user definition only for the Natural/Essential boundary node number
INPUTS ::
NBElementObject = Neuman Boundary nodes object containing 1D element shape functions
fGlobal = Global force vector (double)
NODE_COORD = Node coordinates of the whole domain
Econnectivity = Element connectivities (obtained from readMesh function)
tForce = Traction forces values defined by user (an array of form {Fx, Fy})
NBnodes = Containing node connectivities of elements on the neuman boundary
t = Thickness of the domain
mesh = mesh object containg mesh details and neuman boundary edge nodes
*/
void ImposeTraction(ShapeFn1D NBElementObject, Eigen::VectorXd &fGlobal, std::vector<std::vector<double>> NODE_COORD,
                    std::vector<std::vector<unsigned int>> NBnodes, readMesh mesh, Parameters_LE input);