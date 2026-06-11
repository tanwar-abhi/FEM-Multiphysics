#pragma once

#include "MaterialBase.hpp"

// Data structure to hold material properties of domain for Linear Elastic problems.
class LinearElasticMaterial : private MaterialBase
{
    public:
    // Default Constructor and Destructors
    LinearElasticMaterial();
    ~LinearElasticMaterial();

    LinearElasticMaterial(const float& youngModulus, const float& poissonRatio, 
                            const double& angularVelocity=0, const double& density
    );

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

    // // Static function so that it can be used without an object, this function reads all material properties into a vector,
    // each element is a different material, defined by user in material input file
    std::vector<LinearElasticMaterial> readMaterialInputs();

    Eigen::MatrixXd HookesLawTensor(int eqnType);
};
