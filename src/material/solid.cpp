
#include "solid.hpp"

Solid::Solid() {

}

Solid::~Solid() {

}

Solid::Solid(const float& E, const float& NU, const double& RHO, const double& angularVelocity=0)
{
    youngsModulus = E;
    poissonRatio = NU;
    materialDensity = RHO;
    if (angularVelocity != 0)
    {
        omega = angularVelocity;
    }
}



const std::unordered_map<std::string, MaterialType> materialMap = {
    {"LINEARELASTIC", MaterialType::LinearElastic},
    {"HYPERELASTIC", MaterialType::Hyperelastic},
    {"NONLINEAR", MaterialType::NonLinear},
    {"ELASTOPLASTIC", MaterialType::ElastoPlastic},
    {"VISCOPLASTICITY", MaterialType::Viscoplasticity},
    {"UNKNOWN", MaterialType::Unknown}
};



MaterialType Solid::lookUpMaterialType(const std::string& materialFromJson)
{
    std::transform(materialFromJson.begin(), materialFromJson.end(), materialFromJson.begin(), ::toupper);
    auto it = materialMap.find(materialFromJson);
    if (it != materialMap.end()){
        return it->second;
    }
    std::cout<<"[debug] Solid::lookUpMaterialType material type is not found in lookup table.\n";
    return MaterialType::Unknown;
}

MaterialType Solid::getMaterialType()
{
    return _materialType;
}


// Calculate constitutive equation Hookes Law Tensor for equation
Eigen::MatrixXd Solid::HookesLawTensor(const EquationType& equationType)
{
    // Initializing Linear Elasticity matrix
    Eigen::MatrixXd DMatrix;

    switch (equationType)
    {
        case EquationType::PlaneStress:
        {
            DMatrix = Eigen::MatrixXd::Zero(3, 3);
            DMatrix.row(0) << 1, poissonRatio, 0;
            DMatrix.row(1) << poissonRatio, 1, 0;
            DMatrix.row(2) << 0, 0, (1-poissonRatio)/2;

            double term = youngsModulus/(1 - (poissonRatio*poissonRatio));
            DMatrix *= term;

            break;
        }
        case EquationType::PlaneStrain:
        {
            DMatrix = Eigen::MatrixXd::Zero(3, 3);
            DMatrix.row(0) << (1-poissonRatio), poissonRatio, 0;
            DMatrix.row(1) << poissonRatio, (1-poissonRatio) , 0;
            DMatrix.row(2) << 0, 0, (1-2*poissonRatio)/2;

            double term = youngsModulus / ((1+poissonRatio)*(1-2*poissonRatio));
            DMatrix *= term;

            break;
        }
        case EquationType::Axisymmetric2D:
        {
            DMatrix = Eigen::MatrixXd::Zero(4, 4);
            DMatrix.row(0) << 1-poissonRatio , poissonRatio , poissonRatio , 0;
            DMatrix.row(1) << poissonRatio , 1-poissonRatio , poissonRatio , 0;
            DMatrix.row(2) << poissonRatio , poissonRatio , 1-poissonRatio , 0;
            DMatrix.row(3) << 0 , 0 , 0 , (1-2*poissonRatio)/2.0;

            double term = youngsModulus / ((1+poissonRatio)*(1-2*poissonRatio));
            DMatrix *= term;

            break;
        }
        case EquationType::LinearElastic3D:
        {
            DMatrix = Eigen::MatrixXd::Zero(6, 6);
            DMatrix.row(0) << 1, poissonRatio/(1-poissonRatio), poissonRatio/(1-poissonRatio), 0, 0, 0;
            DMatrix.row(1) << poissonRatio/(1-poissonRatio), 1, poissonRatio/(1-poissonRatio), 0, 0, 0;
            DMatrix.row(2) << poissonRatio/(1-poissonRatio), poissonRatio/(1-poissonRatio), 1, 0, 0, 0;
            DMatrix.row(3) << 0, 0, 0, (1-2*poissonRatio)/(2*(1-poissonRatio)), 0, 0;
            DMatrix.row(4) << 0, 0, 0, 0, (1-2*poissonRatio)/(2*(1-poissonRatio)), 0;
            DMatrix.row(5) << 0, 0, 0, 0, 0, (1-2*poissonRatio)/(2*(1-poissonRatio));

            double term = youngsModulus*(1-poissonRatio)/((1+poissonRatio)*(1-2*poissonRatio));
            DMatrix *= term;

            break;
        }
        
        default:
        {
            std::cerr<<"Material Error : Invalid linear ealstic Material properties\n";
            break;
        }
    }

    return DMatrix;
}