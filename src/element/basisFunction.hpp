
#pragma once

#include "quadratures.hpp"

// #include "../Libraries/Eigen/Core"

// #include <string>
// #include <cmath>
#include <Eigen/Core>

// The internal structure for Basis Function creation, used internally by the main classes i.e. primarily 
// used by ShapeFunction1D and ShapeFunction2d, this only call upon respective 1d or 2d shape function classes
struct BasisFunction{

    BasisFunction();
    ~BasisFunction();
    BasisFunction(int, int);

    // No. of gauss points
    int NGP;

    // elemenType -> Based on the element type definition in gmsh file format
    // Check this URL for more information -> https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    int Etype;

    // Weigths of gauss points
    std::vector<double> wgp;

    // Shape functions for 1D element
    // ngp = no. og gauss points.
    // Returns a vector containing Shape functions for 1D element, return vecctor format ::
    // returnVector[0] = N; returnVector[1] = Nxi; returnVector[2] = N2xi; 
    std::vector<Eigen::MatrixXd> ShapeFunction_1D(const int ngp);

    // Shape function for 2D element
    // Returns a vector containing Shape functions for 2D element
    std::vector<Eigen::MatrixXd> ShapeFunction_2D(const int NGP);

};


class ShapeFunction1D
{
    public:
    int eType;
    std::vector<double> WGP;

    ShapeFunction1D();

    // Since we are calling virtual function, thus we need to define virtual distructor
    virtual ~ShapeFunction1D();

    // No. of gauss points
    int NGP;

    /* Function to calculate 1d shape functions and bind them to element of class.
    Inputs :: 
        npg = No. of gauss points ; 
        elementType = element type identifier as per .msh
    */
    // virtual void getShapeFunction(int ngp, int elementType);
    virtual void getShapeFunction(int elementType);

    // Function to get shape functions for specific gauss points (used in post processing and special cases)
    virtual void getShapeFunction(int elementType, int ngp);

    Eigen::MatrixXd N;
    Eigen::MatrixXd Nxi;
    Eigen::MatrixXd N2xi;
};


// 2D shape functions class inheriting from 1D shape function class
class ShapeFunction2D : public ShapeFunction1D
{
    public:
    ShapeFunction2D();
    ~ShapeFunction2D();

    // Generate 2D shape functions and bind them to element of class
    // Inputs required are npg: no. of gauss points and element type as per .msh
    // void getShapeFunction(int ngp, int elementType);
    void getShapeFunction(int elementType);

    // Function to get shape functions for specific gauss points (used in post processing or special cases)
    void getShapeFunction(int elementType, int ngp);

    Eigen::MatrixXd Neta;
    Eigen::MatrixXd N2eta;
    Eigen::MatrixXd N2xieta;
    Eigen::MatrixXd N2etaxi;
};




// 3D shape functions class inheriting from 2D shape function class
class ShapeFunction3D : public ShapeFunction2D
{
    public:
    ShapeFunction3D();
    ~ShapeFunction3D();

    // Differentiation of shape function w.r.t. 3rd coordinate natural coordinate {xi,eta,zeta}
    Eigen::MatrixXd Nzeta; Eigen::MatrixXd N2zeta;
    Eigen::MatrixXd N2ZetaXi; Eigen::MatrixXd N2ZetaEta;
    Eigen::MatrixXd N2XiZeta; Eigen::MatrixXd N2EtaZeta;

    // Function to get shape functions for specific element types with defined gauss points
    void getShapeFunction(int elementType, int ngp);

};
