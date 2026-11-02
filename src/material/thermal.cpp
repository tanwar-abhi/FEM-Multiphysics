
#include "thermal.hpp"


MaterialThermal::MaterialThermal()
{
    std::cout<<"[INFO] intialized thermal material constructor\n";
}


MaterialThermal::~MaterialThermal()
{

}


Eigen::MatrixXd MaterialThermal::getMaterialThermalConductivity(const int& dimension)
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




std::vector<std::optional<MaterialThermal>> MaterialThermal::readMaterialInputs()
{

    // Temporary object to read material properties defined in material input file
    std::optional<MaterialThermal> material;

    // Vector containing multiple material definations
    std::vector<std::optional<MaterialThermal>> vectorMaterialThermal;

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

        vectorMaterialThermal.push_back(material);
    }

    return vectorMaterialThermal;
}
