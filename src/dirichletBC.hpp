
#pragma once

#include "boundaryInputs.hpp"
#include "element.hpp"
#include "solverInputs.hpp"

// #include <map>
#include <Eigen/SparseCore>

#include <set>

typedef Eigen::SparseMatrix<double> SpMatDouble;
typedef Eigen::Triplet<double> Tr;




/* Apply all dirichlet boundary conditions defined by user in boundary.json and generate reduced stiffness matrix. 
This should be called irrespective of dirichlet boundary being present in problem or not, as it will generate global stiffness matrix

Returns:: Sparse Matrix of Global stiffness, i.e. reduced stiffness matrix after elimination
Inputs::
    tripletGlobal = Triplet containg row,column and Value for Global stiffness Matrix,
    fGlobal = Global RHS load vector
    meshElement[] = Array of element class containing all elemenent connectivities and shape function
    boundary = object of Boundary conditions containing associated parameters
    solverInp = Object of SolverInp class containing sovler specific information
    numNodesDomain = Total number of nodes in the domain
*/
Eigen::SparseMatrix<double> applyDirichletBC(std::vector<Eigen::Triplet<double>> &tripletGlobal, Eigen::VectorXd &fGlobal, const Element meshElement[], const BoundaryConditions &boundary, const SolverInp &solverInp, unsigned long numNodesDomain);




/* Function to impose the dirichlet values at the nodes in the solution to get final solution
Returns :: Final solution in form of Eigen::MatrixXd
Inputs::
    solutionRed = Reduced solution vector obtained eliminating the dirichlet nodes
    meshElement[] = Array of element class containing all elemenent connectivities and shape function
    boundary = object of Boundary conditions containing associated parameters
    numNodesDomain = Total number of nodes in the domain
    dof = Degree of freedom at each node of the domain
*/
Eigen::MatrixXd imposeDirichletValues(const Eigen::VectorXd &solutionRed, const Element meshElement[], const BoundaryConditions &boundary, unsigned long numNodesDomain, int dof);

