
#pragma once
#include <vector>
#include <iostream>

struct Quadrature1D{

    // Gauss Point of quadrature
    double zgp1;

    // Gauss weights
    double wgp;

    Quadrature1D();
    virtual ~Quadrature1D();

    // std::vector<Quadrature1D> GaussPoints2D(const int);

    // Input argument is number of gauss points.
    virtual std::vector<Quadrature1D> GaussQuadrature(const int NGAUS);
};


// 2d Gauss Quadratures
struct Quadrature2D : Quadrature1D{

    Quadrature2D();
    ~Quadrature2D();

    // Second Gauss Point of 2D Quad1D.
    double zgp2;

    // ngp: No. of gauss points; elementType of the mesh
    // elementType is as per GMSH documentation : https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    std::vector<Quadrature2D> GaussQuadrature(const int ngp, int elementType);
};


// 3D Gauss Quadratures
struct Quadrature3D : Quadrature2D{

    // Constructors and Destrucrtors
    Quadrature3D();
    ~Quadrature3D();

    // Third Gauss Point for 3D gauss Points
    double zgp3;

    // Function to call Gauss Quadratures
    // NGP : number of gauss points; 
    // elemType : Type of 3D element as per GMSH documentation : https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    std::vector<Quadrature3D> GaussQuadrature(const int NGP, int elementType);
};
