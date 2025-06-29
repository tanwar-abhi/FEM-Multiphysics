// Assembly routine for the Global system of equation

#pragma once
// #include "BasisFunction.hpp"
// #include "readMesh.hpp"
#include "boundary_HT.hpp"
// #include "Material.hpp"
// #include "SolverInputs.hpp"
// #include "BoundaryInputs.hpp"
#include <Eigen/SparseCore>
#include <Eigen/LU>




/* Function to assemble the global matrix for Linear Elasticity Problem
Input Arguments ::
    tripletStiffness = Triplet containg row,column and Value for Global stiffness Matrix,
    tripletMass = Triplet containing row,column and Value for Global mass matrix,
    f_global = Global foces and load vector
    NODE_COORD = Nodal coordinates of domain,
    meshElement = Array of element class containing domain elemenent connectivities and shape function
    domainElemID = Domian element position within meshElement
    material = object of linearElasticMaterial class
    boundary = object of BoundaryConditions that contain all constants and type of variable
    eqnIndex = equation number defined in the solver inputs
*/
void assembleMatrixLE(std::vector<Tr> &tripletStiffness, Eigen::VectorXd &f_global, const std::vector<std::vector<double>> &NODE_COORD, const Element &meshElement,
                    int domainElemID, const SolverInp &solverInp, std::vector<LinearElasticMaterial> material, const BoundaryConditions &boundary, int eqnIndex);





/* Calculate jacobian matrix for 2D element, jacobian is used to trasnform local integral to parametric coordinates
    This function return the "Jacobian matrix" in form of Eigen::Matrix2D i.e. 2x2 matrix
    Nx = Derivative of shape functioin w.r.t. x
    Ny = Derivative of shape function w.r.t. y
    Xe = Coordinate of nodes of current element
    BasisFn = Object of shapeFn containing all shape functions and it's derivatives in natural coordinate(w.r.t xi,eta)
    i = iterator of the current loop of gauss point
*/
Eigen::Matrix2d Jacobian2D(Eigen::VectorXd &Nx, Eigen::VectorXd &Ny, const Eigen::MatrixXd &Element_NC, const ShapeFunction2D &BasisFn, const int i);



/* Calculate jacobian matrix for 3D element, it's used to transform local integral to isotropic coordinates and vice versa
    This function return the "Jacobian matrix" in form of Eigen::Matrix2D i.e. 3x3 matrix
    Nx = Derivative of shape functioin w.r.t. x
    Ny = Derivative of shape function w.r.t. y
    Nz = Derivative of shape function w.r.t. z
    Xe = Coordinate of nodes of current element {Element_NC}
    shapeFn = objecto of 3D shapefunction for the respective element
    i = iterator of the current looop of gauss point
*/
Eigen::Matrix3d Jacobian3D(Eigen::VectorXd &Nx, Eigen::VectorXd &Ny, Eigen::VectorXd &Nz, const Eigen::MatrixXd &Xe, const ShapeFunction3D &shapeFn, const int itr);




/* Function to assemble the global matrix for Heat Transfer Problems
Input Arguments ::
    tripletStiffness = Triplet containg row,column and Value for Global stiffness Matrix,
    tripletMass = Triplet containing row,column and Value for Global mass matrix,
    f_global = Global thermal load vector
    NODE_COORD = Global Nodal coordinates of domain,
    meshElement = Array of element class containing domain elemenent connectivities and shape function
    domainElemId = Element id of the domain elements
    solverInp = Object of SolverInp class containing sovler specific information
    material = Thermal material properties of the material
    boundary = object of Boundary conditions containing associated parameters
*/
void assembleMatrixHT(std::vector<Tr> &tripletStiffness, std::vector<Tr> &tripletMass, Eigen::VectorXd &f_global, const std::vector<std::vector<double>> &NODE_COORD,
                     const Element meshElement[], int domainElemID, const SolverInp &solverInp, const std::vector<MaterialThermal> &material, const BoundaryConditions &boundary);



