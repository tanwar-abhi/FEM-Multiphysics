#pragma once

#include <eigen3/Eigen/Core>

#include "solverInputs.hpp"
#include "MaterialBase.hpp"


enum class MaterialType
{
    LinearElastic,
    Hyperelastic,
    NonLinear,
    ElastoPlastic,
    Viscoplasticity,
    Unknown
};


// Data structure to hold material properties of domain for Linear Elastic problems.
class Solid : public MaterialBase<Solid>
{
    public:
    // Default Constructor and Destructors
    Solid();
    ~Solid();

    // Inherit the constructior from Base class.
    Solid(const std::string& inputsDirectoryPath) : MaterialBase(inputsDirectoryPath){};

    Solid(const float& E, const float& NU, const double& density, const double& angularVelocity=0);

    // // Static function so that it can be used without an object, this function reads all material properties into a vector,
    // each element is a different material, defined by user in material input file
    std::vector<Solid> readMaterialInputs();

    Eigen::MatrixXd HookesLawTensor(const EquationType& eqnType);

    MaterialType getMaterialType();

    // Thickness of the domain (only for 2D problems)
    float thickness = 1.0;

    // Young's modulus [E] in (Mpa)
    float youngsModulus = 0.0;

    // Poisson's Ratio [Nu]
    float poissonRatio = 0.0;

    // Angular velocity {only for Axisymmetric case}
    double omega = 0.0;

    // Material density [RHO] {only for Axisymmetric case} 
    double materialDensity = 0;
    
    private:

    // Type of the material
    MaterialType _materialType = MaterialType::Unknown;

    MaterialType lookUpMaterialType(const std::string& materialString);
};
