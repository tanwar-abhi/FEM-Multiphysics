#include "MaterialBase.hpp"

MaterialBase::MaterialBase()
{

}

MaterialBase::~MaterialBase()
{

}

void MaterialBase::setName(const std::string& value)
{
    name = value;
}

std::string MaterialBase::getName()
{
    return name;
}

void MaterialBase::setElementTagId(const int& tagID)
{
    elementTagId = tagID;
}

int MaterialBase::getElementTagId()
{
    return elementTagId;
}


template <typename T>
std::vector<T> MaterialBase::readMaterialInputs()
{

    // Temporary object to read material properties defined in material input file
    T material;

    // Vector containing multiple material definations
    std::vector<T> MaterialThermal_collection;

    std::ifstream materialText(_directoryInputs + "solver.json");
    Json::Value mat_root;
    Json::Reader mat_reader;
    bool parsingSuccessful = mat_reader.parse(materialText, mat_root);
    if(!parsingSuccessful)
    {
        std::cout << "Error parsing the string" << std::endl;
    }
    Json::Value inp_materialProp = mat_root["materialProperty"];

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

        MaterialThermal_collection.push_back(material);
    }

    return MaterialThermal_collection;
}

