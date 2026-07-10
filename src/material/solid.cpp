
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


// Read material data from JSON input files
std::vector<Solid> Solid::readMaterialInputs()
{
    Solid material;
    std::vector <Solid> materialVector;

    const std::string solverInputJsonFile = getInputsDirectory() + getMaterialInputFileName();
    std::ifstream mat_text(solverInputJsonFile);
    Json::Value mat_root;
    Json::Reader mat_reader;
    bool parsingSuccessful = mat_reader.parse(mat_text, mat_root);

    if(!parsingSuccessful)
    {
	    std::cerr << "[Error] Solid::readMaterialInputs : Error unable to parse the string" << std::endl;
        std::runtime_error("Error unable to parse the string in Solid::readMaterialInputs");
    }

    const Json::Value materialJson = mat_root["materialProperty"];

    // Temporary Variable used for reading inputs
    std::string tmp;

    if (materialJson.isArray())
    {
        for(int index = 0, n = materialJson.size(); index < n; index++)
        {
            if (materialJson[index].isMember("type") && materialJson[index]["type"].isString())
            {
                tmp = materialJson[index]["type"].asString();
                material._materialType = lookUpMaterialType(tmp);
            }
    
            if (materialJson[index].isMember("name") && materialJson[index]["name"].isString())
            {
                tmp = materialJson[index]["name"].asString();
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
                material.setName(tmp);
            }
    
            if (materialJson[index].isMember("rho") && materialJson[index]["rho"].isDouble())
            {
                material.materialDensity = materialJson[index]["rho"].asDouble();
            }

            if (materialJson[index].isMember("nu") && materialJson[index]["nu"].isDouble())
            {
                material.poissonRatio = materialJson[index]["nu"].asDouble();
            }

            if (materialJson[index].isMember("thickness") && materialJson[index]["thickness"].isDouble())
            {
                material.thickness = materialJson[index]["thickness"].asDouble();
            }

            if (materialJson[index].isMember("youngsModulus") && materialJson[index]["youngsModulus"].isDouble())
            {
                material.youngsModulus = materialJson[index]["youngsModulus"].asDouble();
            }

            if (materialJson[index].isMember("angularVelocity") && materialJson[index]["angularVelocity"].isDouble())
            {
                material.omega = materialJson[index]["angularVelocity"].asDouble();
            }
            else{
                if (material.getMaterialType() != MaterialType::LinearElastic)
                {
                    std::cerr<<"[ERROR] Solid::readMaterialInputs";
                    std::runtime_error("Material type is unavailable, currently only support LinearElastic materials");
                }
            }

            // Eliminate user input error for thickness
            if (material.thickness <= 0)
            {
                material.thickness = 1;
            }
    
            materialVector.emplace_back(material);
        }
    }
    else{
        std::runtime_error("[ERROR] Solid::readMaterialInputs solverConfig.json file does't contain materialProperty array." );
    }

    return materialVector;
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
            double term = youngsModulus/(1 - pow(poissonRatio,2));
            DMatrix.row(0) << 1, poissonRatio, 0;
            DMatrix.row(1) << poissonRatio, 1, 0;
            DMatrix.row(2) << 0, 0, (1-poissonRatio)/2;
            DMatrix *= term;
            break;
        }
        case EquationType::PlaneStrain:
        {
            DMatrix = Eigen::MatrixXd::Zero(3, 3);
            double term = youngsModulus / ((1+poissonRatio)*(1-2*poissonRatio));
            DMatrix.row(0) << (1-poissonRatio), poissonRatio, 0;
            DMatrix.row(1) << poissonRatio, (1-poissonRatio) , 0;
            DMatrix.row(2) << 0, 0, (1-2*poissonRatio)/2;
            DMatrix *= term;
            break;
        }
        case EquationType::Axisymmetric2D:
        {
            DMatrix = Eigen::MatrixXd::Zero(4, 4);
            double term = youngsModulus / ((1+poissonRatio)*(1-2*poissonRatio));
            DMatrix.row(0) << 1-poissonRatio , poissonRatio , poissonRatio , 0;
            DMatrix.row(1) << poissonRatio , 1-poissonRatio , poissonRatio , 0;
            DMatrix.row(2) << poissonRatio , poissonRatio , 1-poissonRatio , 0;
            DMatrix.row(3) << 0 , 0 , 0 , (1-2*poissonRatio)/2.0;
            DMatrix *= term;
            break;
        }
        case EquationType::LinearElastic3D:
        {
            DMatrix = Eigen::MatrixXd::Zero(6, 6);
            double term = youngsModulus*(1-poissonRatio)/((1+poissonRatio)*(1-2*poissonRatio));
            DMatrix.row(0) << 1, poissonRatio/(1-poissonRatio), poissonRatio/(1-poissonRatio), 0, 0, 0;
            DMatrix.row(1) << poissonRatio/(1-poissonRatio), 1, poissonRatio/(1-poissonRatio), 0, 0, 0;
            DMatrix.row(2) << poissonRatio/(1-poissonRatio), poissonRatio/(1-poissonRatio), 1, 0, 0, 0;
            DMatrix.row(3) << 0, 0, 0, (1-2*poissonRatio)/(2*(1-poissonRatio)), 0, 0;
            DMatrix.row(4) << 0, 0, 0, 0, (1-2*poissonRatio)/(2*(1-poissonRatio)), 0;
            DMatrix.row(5) << 0, 0, 0, 0, 0, (1-2*poissonRatio)/(2*(1-poissonRatio));
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