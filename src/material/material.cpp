
#include "material.hpp"

Material::Material()
{
	// Default constructor of Material class
}

Material::~Material()
{
	// Default Destructor
}



bool Material::hasSolid() const
{
	return solid.has_value();
}

bool Material::hasThermal() const
{
	return thermal.has_value();
}




/*
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

            if (material.getMaterialType() != MaterialType::LinearElastic)
            {
                std::cerr<<"[ERROR] Solid::readMaterialInputs";
                std::runtime_error("Material type is unavailable, currently only support LinearElastic materials");
            }

            // Eliminate user input error for thickness, or missing information
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
*/





/*
std::vector<MaterialThermal> MaterialThermal::readMaterialInputs()
{

    // Temporary object to read material properties defined in material input file
    MaterialThermal material;

    // Vector containing multiple material definations
    std::vector<MaterialThermal> vectorMaterialThermal;

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
*/