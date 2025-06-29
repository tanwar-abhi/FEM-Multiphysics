
#pragma once
// #include "ProblemParameters.hpp"
#include <Boundary_HT.hpp>
#include <DirichletBC.hpp>
#include <PostProcess.hpp>

// Time discretization function using theta method for heat transfer problems
// void TimeDiscretization_ThetaHT(Eigen::MatrixXd &Temperature ,SpMatDouble C, SpMatDouble K, Eigen::VectorXd fglobal,std::vector<std::vector<double>> NODE_COORD, ProbParameters_HT inputs, std::vector<std::vector<unsigned int>> T0nodes, readMesh mesh, std::set<unsigned int> DBnodes);
void TimeDiscretization_ThetaHT(Eigen::MatrixXd &Temperature ,SpMatDouble C, SpMatDouble K, Eigen::VectorXd fglobal, std::vector<std::vector<unsigned int>> Econnectivity,std::vector<std::vector<double>> NODE_COORD, std::vector<std::vector<unsigned int>> T0nodes, readMesh mesh, std::set<unsigned int> DBnodes, SolverInp solverInp, NeumannBC boundary, MaterialThermal mat);



/* Impose initial boundary conditions for the first time step
temp -> Nodal Temperatures vector 
nodeBoundary -> Node numbers at the boundary where initial condition must be imposed
value -> Initial Temperature value to be imposed
*/
void initialConditions(Eigen::VectorXd &temp, std::vector<unsigned int> nodeBoundary, double value);
