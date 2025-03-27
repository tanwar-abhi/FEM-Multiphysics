
#pragma once
#include <vector>
#include <string>
#include <Eigen/Core>
#include <algorithm>
#include <fstream>
#include <json/value.h>
#include <json/json.h>
#include <iostream>


// Data structure to hold material properties of domain for Linear Elastic problems.
class LinearElasticMaterial
{
    public:
    // Default Constructor and Destructors
    LinearElasticMaterial();
    ~LinearElasticMaterial();

    // Material Defination identifier name
    std::string name;

    // Thickness of the domain (only for 2D problems)
    float thickness = 1.0;

    // Young's modulus in (Mpa)
    float E = 0;

    // Poisson's Ratio
    float NU = 0;

    // Angular velocity {only for Axisymmetric case}
    double omega = 0;

    // Material density {only for Axisymmetric case}
    double RHO = 0;

    // Element tag id to map material with equation and boundary
    int elemTagId;

    // Static function so that it can be used without an object, this function reads all material properties into a vector,
    // each element is a different material, defined by user in materil input file
    static std::vector<LinearElasticMaterial> readMaterialInputs();

    Eigen::MatrixXd HookesLawTensor(int eqnType);
};


// Data structure to hold material properties of domain for heat transfer problem.
class MaterialThermal{

    public:
    // Default Constructor and Destructors
    MaterialThermal();
    ~MaterialThermal();

    // Thickness of domain (element) to be defined for 2D cartesian coordinate problem only {units : m}
    // NOTE -> Never set thickness = 0, by default value of thickenss = 1
    double thickness = 1;

    // Material Defination identifier name
    std::string name;

    // Mass Density {SI units = kg/m3}
    double RHO = 0;

    // Specific heat (c) {SI units in J/Kg °C}
    double spHeat = 0;

    /* Thermal conductivity for isotropic material {SI unit = W/(m °C)}
    k would be a tensor value for anisotropic material 
    In 2D, k = [Kxx, Kxy; Kyx, Kyy]; 
    In 3D k = [Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz]
    */
    double k[4] = {0, 0, 0, 0};

    // Element tag id to map material with equation and boundary
    int elemTagId;

    // Static function can be used without object, reads all material properties for Heat Transfer defined in inputs file
    static std::vector<MaterialThermal> readMaterialInputs();
};




// Returns a matrix of thermal conductivity for anisotropic material
Eigen::MatrixXd Material_HT(int dimension,const double *k);


// Returns a matrix of thermal conductivity for isotropic material
// Eigen::MatrixXd Material_HT(int dimension, double K);
