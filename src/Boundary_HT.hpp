// The functions defined within this function are to impose contribution of the heat transfer terms
// for the heat transfer problems.

#pragma once
#include "BasisFunction.hpp"
#include "readMesh.hpp"
#include "SolverInputs.hpp"
#include "BoundaryInputs.hpp"
#include "Material.hpp"
#include "Element.hpp"
#include <Eigen/SparseCore>
#include <Eigen/Geometry>


typedef Eigen::SparseMatrix<double> SpMatDouble;
typedef Eigen::Triplet<double> Tr;



// Heat load vector contribution due to line source/sink
void LineHeatSource(Eigen::VectorXd &fGlobal, std::vector<std::vector<double>> NODE_COORD, readMesh mesh, std::vector<std::vector<unsigned int>> T, SolverInp solverInp, MaterialThermal object, NeumannBC input);


/*
Inputs::
BasisFn1D : Object of ShapeFn1D class containig all 1D shape functions and their derivatives
fGlobal : Gloabal force vector
NODE_COORD : Nodal points coordinates of domain
mesh : Object of readMesh class containg information about the mesh
input : Object of ProbParameter_HT class containg all the given input by the user
RadBNodes : a 2D vector containing connectivities of all nodes at the radiant heat flow boundary

Heat load vector contribution due to radiant heat flow */
// void SurfaceRadiation_HF(ShapeFn1D BasisFn1D, Eigen::VectorXd &fGlobal, std::vector<std::vector<double>> NODE_COORD, Element meshElement, std::vector<std::vector<unsigned int>> RadBNodes, SolverInp SIobject, NeumannBC boundary, MaterialThermal mat);




/* Calculates the unit outward normal from the surface of element
Returns a 2D vector containg unit vectors in form [nx, ny] for 2D problems and [nx,ny,nz] for 3D problems

INPUTS::
nodeCoord = Element node coordintes
elemLength = Length of element
elemType = Type of element
*/
std::vector<double> outwardNormal(const Eigen::MatrixXd &nodeCoord, double elemLength, int elemType);





// ################################## NEW improved #############
/*
Conduction boundary contributions to global thermal load vector.
Inputs ::
fGlobal : Gloabal force vector
NODE_COORD : Nodal points coordinates of domain
meshElement : Element object of elements array, containing Element connectivity (conduction boundary), shape Functions, element type, order and other element details
solverInp : Object of SolverInp class containing data regarding solver type, coordinate system and other details given by user
boundary : object of the Neuman class containing user defined details of condution boundary condition
*/
void SurfaceConduction(Eigen::VectorXd &fGlobal, const std::vector<std::vector<double>> &NODE_COORD, const Element &meshElement, 
                        const SolverInp &solverInp, const NeumannBC &boundary, const MaterialThermal &material);


/*
Convection boundary contribution to global thermal load vector.
Inputs
fGlobal : Gloabal force vector
NODE_COORD : Nodal points coordinates of domain
meshElement : Element object of elements array, containing Element connectivity (convection boundary), shape Functions, element type, order and other element details
solverInp : Object of SolverInp class containing data regarding solver type, coordinate system and other details given by user
boundary : object of the Neuman class containing user defined details of convection boundary condition
*/
void SurfaceConvection(Eigen::VectorXd &fGlobal, const std::vector<std::vector<double>> &NODE_COORD, const Element &meshElement, 
                        const SolverInp &solverInp, const NeumannBC &boundary, const MaterialThermal &material);


/*
Stiffness matrix contribution from convective terms due to convection boundary.
Inputs ::
TrList : Triplet list containing non zero values, row, column number for global stiffness matrix
NODE_COORD : Nodal points coordinates of domain
meshElement : Element object of elements array, containing Element connectivity (convection boundary), shape Functions, element type, order and other element details
solverInp : Object of SolverInp class containing data regarding solver type, coordinate system and other details given by user
boundary : object of the Neuman class containing user defined details of convection boundary condition
*/
void StiffnessConvection(std::vector<Eigen::Triplet<double>> &trList, const std::vector<std::vector<double>> &NODE_COORD,
                            const Element &, const SolverInp &, const NeumannBC &,  const MaterialThermal &material);
