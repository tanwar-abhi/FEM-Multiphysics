
#pragma once
#include <vector>
#include <iostream>

struct Quad1D{

    // Gauss Point of quadrature
    double zgp1;
    
    // Gauss weights
    double wgp;

    Quad1D();
    virtual ~Quad1D();

    // std::vector<Quad1D> GaussPoints2D(const int);

    // Input argument is number of gauss points.
    virtual std::vector<Quad1D> GaussQuad(const int NGAUS);
};


// 2d Gauss Quadratures
struct Quad2d : Quad1D{

    Quad2d();
    ~Quad2d();

    // Second Gauss Point of 2D Quad1D.
    double zgp2;

    // ngp: No. of gauss points; elementType of the mesh
    // elementType is as per GMSH documentation : https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    std::vector<Quad2d> GaussQuad(const int ngp, int elementType);
};


// 3D Gauss Quadratures
struct Quad3D : Quad2d{

    // Constructors and Destrucrtors
    Quad3D();
    ~Quad3D();

    // Third Gauss Point for 3D gauss Points
    double zgp3;

    // Function to call Gauss Quadratures
    // NGP : number of gauss points; 
    // elemType : Type of 3D element as per GMSH documentation : https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    std::vector<Quad3D> GaussQuad(const int NGP, int elementType);
};
