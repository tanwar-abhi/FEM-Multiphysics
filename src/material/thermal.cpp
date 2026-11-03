
#include "thermal.hpp"


Thermal::Thermal()
{
    spdlog::debug("Thermal::Thermal intialized thermal material constructor");
}


Thermal::~Thermal()
{
    spdlog::debug("Thermal::~Thermal destructor called");
}


Eigen::MatrixXd Thermal::getThermalConductivity(const int& dimension)
{
    Eigen::MatrixXd thermalConductivityMatrix(dimension, dimension);

    if (dimension == 2){
        int nonZero = 0;
        for (int i = 0; i < 4; i++){
            if (thermalConductivity2D[i] != 0){
                nonZero++;
            }
        }

        // Isothermal material (one value of thermal conductivity)
        if (nonZero == 1){
            thermalConductivityMatrix.row(0) << thermalConductivity2D[0], 0;
            thermalConductivityMatrix.row(1) << 0, thermalConductivity2D[0];
        }
        // Anisothermal material
        else{
            thermalConductivityMatrix.row(0) << thermalConductivity2D[0], thermalConductivity2D[1];
            thermalConductivityMatrix.row(1) << thermalConductivity2D[2], thermalConductivity2D[3];
        }
    }
    else if (dimension == 3){
        // Isothermal material (one value of thermal conductivity)
        thermalConductivityMatrix.row(0) << thermalConductivity2D[0] , 0 , 0;
        thermalConductivityMatrix.row(1) << 0 , thermalConductivity2D[0] , 0;
        thermalConductivityMatrix.row(2) << 0 , 0 , thermalConductivity2D[0];
    }

    return thermalConductivityMatrix;
}

const std::unordered_map<std::string, ThermalMaterialType> thermalMaterialMap = 
{
    {"INSULATOR", ThermalMaterialType::Insulator},
    {"CONDUCTOR", ThermalMaterialType::Conductor},
    {"THERMAL", ThermalMaterialType::Thermal},
    {"ISOTROPIC", ThermalMaterialType::Isotropic},
    {"ANISOTROPIC", ThermalMaterialType::Anisotropic},
    {"UNKNONW", ThermalMaterialType::Unknown}
};


void Thermal::lookUpMaterialType(const std::string& materialTypeJson)
{
    auto it = thermalMaterialMap.find(materialTypeJson);
    if (it != thermalMaterialMap.end())
    {
        _materialType = it->second;
        spdlog::info("Thermal material lookup resulted in material type = {}", it->first);
    }
    else {
        spdlog::critical("Thermal::lookUpMaterialType counld't find defined material type in the lookup table.");
        std::runtime_error("Material type = " + materialTypeJson + " is not found in the lookup table.");
    }
}


void Thermal::readThermalMaterialInputs(const Json::Value& materialJson, std::vector<Thermal>& vectorThermal)
{

    // Temporary object to read material properties defined in material input file
    Thermal material;

    //Variable used for cross checking if material type and equation type matches;
    std::string tmp;

    tmp = materialJson["type"].asString();
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
    material.lookUpMaterialType(tmp);

    // Thermal conductivity is of form [Kxx, Kxy, Kyx, Kyy] for anisotropic material
    if (material._materialType == ThermalMaterialType::Anisotropic)
    {
        for (int i = 0; i < 4; i++)
        {
            material.thermalConductivity2D[i] = materialJson["thermalConductivity"][i].asDouble();
        }
    }
    else
    {
        material.thermalConductivity3D[0] = materialJson["thermalConductivity"].asDouble();
    }

    if (materialJson.isMember("name") && materialJson["name"].isString())
    {
        tmp = materialJson["name"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        material._name = tmp;
    }

    if (materialJson.isMember("rho") && materialJson["rho"].isDouble())
    {
        material.massDensity = materialJson["rho"].asDouble();
    }

    if (materialJson.isMember("specificHeat") && materialJson["specificHeat"].isDouble())
    {
        material.specificHeat = materialJson["specificHeat"].asDouble();
    }

    material.thickness = materialJson["thickness"].asDouble();

    // Eliminate user input error for thickness
    if (material.thickness <= 0){
        material.thickness = 1;
        spdlog::info("Thermal::readThermalMaterialInputs detected incorrect material thickness {}", material.thickness);
    }

    vectorThermal.emplace_back(material);
}



/*
void Thermal::readThermalMaterialInputs()
{

    // Temporary object to read material properties defined in material input file
    std::optional<Thermal> material;

    // Vector containing multiple material definations
    std::vector<std::optional<Thermal>> vectorThermal;

    std::ifstream materialJsonFile(getInputsDirectory() + getMaterialInputFileName());
    Json::Value materialRoot;
    Json::Reader materialReader;
    bool parsingSuccessful = materialReader.parse(materialJsonFile, materialRoot);
    if(!parsingSuccessful)
    {
        std::cout << "Error parsing the string" << std::endl;
    }
    Json::Value inp_materialProp = materialRoot["materialProperty"];

    //Variable used for crosschecking if material type and equation type matches;
    std::string tmp;

    for(int index = 0; index < inp_materialProp.size(); index++){

        tmp = inp_materialProp[index]["type"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

        // Thermal conductivity is of form [Kxx, Kxy, Kyx, Kyy] for anisotropic material
        if (tmp == "HEATTRANSFER-ANISO"){
            for (int i = 0; i < 4; i++){
                material.k[i] = inp_materialProp[index]["thermalConductivity"][i].asDouble();
            }
        }
        else{
            material.k[0] = inp_materialProp[index]["thermalConductivity"].asDouble();
        }

        tmp = inp_materialProp[index]["name"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        material.name = tmp;

        material.RHO = inp_materialProp[index]["rho"].asDouble();
        material.spHeat = inp_materialProp[index]["specificHeat"].asDouble();
        material.thickness = inp_materialProp[index]["thickness"].asDouble();

        // Eliminate user input error for thickness
        if (material.thickness <= 0){
            material.thickness = 1;
        }

        vectorThermal.push_back(material);
    }

    return vectorThermal;
}
*/