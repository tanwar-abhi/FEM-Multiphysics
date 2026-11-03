
#include "solid.hpp"

Solid::Solid() {

}

Solid::~Solid() {

}

// Solid::Solid(const float& E, const float& NU, const double& RHO, const double& angularVelocity=0)
// {
//     youngsModulus = E;
//     poissonRatio = NU;
//     materialDensity = RHO;
//     if (angularVelocity != 0)
//     {
//         omega = angularVelocity;
//     }
// }



const std::unordered_map<std::string, SolidMaterialType> materialMap = {
    {"LINEARELASTIC", SolidMaterialType::LinearElastic},
    {"HYPERELASTIC", SolidMaterialType::Hyperelastic},
    {"NONLINEAR", SolidMaterialType::NonLinear},
    {"ELASTOPLASTIC", SolidMaterialType::ElastoPlastic},
    {"VISCOPLASTICITY", SolidMaterialType::ViscoElastic},
    {"VISCOELASTIC", SolidMaterialType::ViscoPlastic},
    {"UNKNOWN", SolidMaterialType::Unknown}
};



void Solid::lookUpMaterialType(const std::string& materialFromJson)
{
    auto it = materialMap.find(materialFromJson);
    if (it != materialMap.end()){
        _materialType = it->second;
        spdlog::info("Solid material lookup resulted in material type = {}", it->first);
    }
    else {
        spdlog::critical("Solid::lookUpMaterialType material type in json file={} , not found in lookup table.", materialFromJson);
        std::runtime_error("material type=" + materialFromJson + " is not found in lookup table.");
    }
}

SolidMaterialType Solid::getMaterialType()
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


void Solid::readSolidMaterialInputs(const Json::Value& materialJson, std::vector<Solid> &materialVector)
{
    std::string tmp;
    Solid material;

    if (materialJson.isMember("type") && materialJson["type"].isString())
    {
        tmp = materialJson["type"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        material.lookUpMaterialType(tmp);
    }

    // if (materialJson.isMember("name") && materialJson["name"].isString())
    // {
    //     tmp = materialJson["name"].asString();
    //     std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
    //     material->setName(tmp);
    // }

    if (materialJson.isMember("rho") && materialJson["rho"].isDouble())
    {
        material.materialDensity = materialJson["rho"].asDouble();
    }

    if (materialJson.isMember("nu") && materialJson["nu"].isDouble())
    {
        material.poissonRatio = materialJson["nu"].asDouble();
    }

    if (materialJson.isMember("thickness") && materialJson["thickness"].isDouble())
    {
        material.thickness = materialJson["thickness"].asDouble();
    }

    if (materialJson.isMember("youngsModulus") && materialJson["youngsModulus"].isDouble())
    {
        material.youngsModulus = materialJson["youngsModulus"].asDouble();
    }
    else{
        spdlog::critical("Missing solid material stiffness value of youngsModulus");
        std::runtime_error("[Error] missing value of youngs Modulus");
    }

    if (materialJson.isMember("angularVelocity") && materialJson["angularVelocity"].isDouble())
    {
        material.omega = materialJson["angularVelocity"].asDouble();
    }

    if (material.getMaterialType() != SolidMaterialType::LinearElastic)
    {
        std::cerr << "[ERROR] Solid::readMaterialInputs";
        spdlog::error("Material type={} is unavailable, currently only support Linear Elastic.", materialJson["type"].asString());
        std::runtime_error("Material type is unavailable, currently only support LinearElastic materials");
    }

    // Eliminate user input error for thickness, or missing information
    if (material.thickness <= 0)
    {
        material.thickness = 1;
    }

    materialVector.emplace_back(material);
    spdlog::debug("Solid::readSolidMaterialInputs after emplace_back solid materialVector.size = {}", materialVector.size());
}
