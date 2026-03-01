
#include "material.hpp"
#include <cmath>

// Definition of the Material class constructor and destructor
MaterialThermal::MaterialThermal(){

}

MaterialThermal::~MaterialThermal(){

}

LinearElasticMaterial::LinearElasticMaterial() {

}

LinearElasticMaterial::~LinearElasticMaterial() {

}




// Read material data from JSON input files
std::vector<LinearElasticMaterial> LinearElasticMaterial::readMaterialInputs()
{
    LinearElasticMaterial material;
    std::vector <LinearElasticMaterial> material_collection;

    std::ifstream mat_text("solver.json");
    Json::Value mat_root;
    Json::Reader mat_reader;
    bool parsingSuccessful = mat_reader.parse(mat_text, mat_root);

    if(!parsingSuccessful)
    {
	    std::cout << "Material Inputs Error : Error parsing the string" << std::endl;
        exit(-404);
    }

    const Json::Value inp_materialProp = mat_root["materialProperty"];

    // Temporary Variable used for reading inputs
    std::string tmp;

    for(int index = 0; index < inp_materialProp.size(); index++)
    {
        tmp = inp_materialProp[index]["type"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

        tmp = inp_materialProp[index]["name"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        material.name = tmp;

        material.RHO = inp_materialProp[index]["rho"].asDouble();
        material.NU = inp_materialProp[index]["nu"].asFloat();
        material.thickness = inp_materialProp[index]["thickness"].asFloat();
        material.E = inp_materialProp[index]["youngsMod"].asFloat();
        material.omega = inp_materialProp[index]["omega"].asDouble();

        // Eliminate user input error for thickness
        if (material.thickness <= 0)
        {
            material.thickness = 1;
        }

        material_collection.push_back(material);
    }

    return material_collection;
}









// Calculate constitutive equation Hookes Law Tensor for equation
Eigen::MatrixXd LinearElasticMaterial::HookesLawTensor(int equationType)
{
    // Initializing Linear Elasticity matrix
    Eigen::MatrixXd DMatrix;

    // Plane Stress
    if (equationType == 1)
    {
        DMatrix = Eigen::MatrixXd::Zero(3, 3);
        double term = E/(1 - pow(NU,2));
        DMatrix.row(0) << 1, NU, 0;
        DMatrix.row(1) << NU, 1, 0;
        DMatrix.row(2) << 0, 0, (1-NU)/2;
        DMatrix *= term;
    }
    // Plane Strain
    else if (equationType == 13)
    {
        DMatrix = Eigen::MatrixXd::Zero(3, 3);
        double term = E / ((1+NU)*(1-2*NU));
        DMatrix.row(0) << (1-NU), NU, 0;
        DMatrix.row(1) << NU, (1-NU) , 0;
        DMatrix.row(2) << 0, 0, (1-2*NU)/2;
        DMatrix *= term;
    }
    // Axisymmetric
    else if (equationType == -1)
    {
        DMatrix = Eigen::MatrixXd::Zero(4, 4);
        double term = E / ((1+NU)*(1-2*NU));
        DMatrix.row(0) << 1-NU , NU , NU , 0;
        DMatrix.row(1) << NU , 1-NU , NU , 0;
        DMatrix.row(2) << NU , NU , 1-NU , 0;
        DMatrix.row(3) << 0 , 0 , 0 , (1-2*NU)/2.0;
        DMatrix *= term;
    }
    // 3D Linear elastic
    else if (equationType == 14)
    {
        DMatrix = Eigen::MatrixXd::Zero(6, 6);
        double term = E*(1-NU)/((1+NU)*(1-2*NU));
        DMatrix.row(0) << 1, NU/(1-NU), NU/(1-NU), 0, 0, 0;
        DMatrix.row(1) << NU/(1-NU), 1, NU/(1-NU), 0, 0, 0;
        DMatrix.row(2) << NU/(1-NU), NU/(1-NU), 1, 0, 0, 0;
        DMatrix.row(3) << 0, 0, 0, (1-2*NU)/(2*(1-NU)), 0, 0;
        DMatrix.row(4) << 0, 0, 0, 0, (1-2*NU)/(2*(1-NU)), 0;
        DMatrix.row(5) << 0, 0, 0, 0, 0, (1-2*NU)/(2*(1-NU));
        DMatrix *= term;
    }
    else{
        std::cerr<<"Material Error : Invalid linear ealstic Material properties\n";
        exit(-101);
    }

    return DMatrix;
}









// Returns Material's thermal conductivity tensor for Isotropic or anisotropic material
Eigen::MatrixXd Material_HT(int dim, const double *k)
{
    // Thermal conductivity tensor of form [Kxx, Kxy; Kyx, Kyy]
    Eigen::MatrixXd thermalConductivity(dim, dim);

    if (dim == 2){
        int nonZero = 0;
        for (int i = 0; i < 4; i++){
            if (k[i] != 0){
                nonZero++;
            }
        }

        // Isothermal material (one value of thermal conductivity)
        if (nonZero == 1){
            thermalConductivity.row(0) << k[0], 0;
            thermalConductivity.row(1) << 0, k[0];
        }
        // Anisothermal material
        else{
            thermalConductivity.row(0) << k[0], k[1];
            thermalConductivity.row(1) << k[2], k[3];
        }
    }
    else if (dim == 3){
        // Isothermal material (one value of thermal conductivity)
        thermalConductivity.row(0) << k[0] , 0 , 0;
        thermalConductivity.row(1) << 0 , k[0] , 0;
        thermalConductivity.row(2) << 0 , 0 , k[0];
    }

    return thermalConductivity;
}



std::vector<MaterialThermal> MaterialThermal::readMaterialInputs()
{

    // Temporary object to read material properties defined in material input file
    MaterialThermal material;

    // Vector containing multiple material definations
    std::vector<MaterialThermal> MaterialThermal_collection;

    std::ifstream mat_text(_directoryInputs + "solver.json");
    Json::Value mat_root;
    Json::Reader mat_reader;
    bool parsingSuccessful = mat_reader.parse(mat_text, mat_root);
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



void MaterialThermal::setInputDirectory(const std::string& directoryPath)
{
    _directoryInputs = directoryPath;
}


std::string MaterialThermal::getInputDirectory()
{
    return _directoryInputs;
}

