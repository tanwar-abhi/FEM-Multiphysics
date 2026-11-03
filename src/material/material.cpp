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
    return !solid.empty();
}

bool Material::hasThermal() const
{
	return !thermal.empty();
}

std::string Material::getMaterialInputFileName() const
{
    return _materialInputFileName;
}

void Material::setInputsDirectory(const std::string& inputsFolderPath)
{
    _directoryInputs = inputsFolderPath;
}

std::string Material::getInputsDirectory() const
{
    return _directoryInputs;
}

void Material::setName(const std::string& newName)
{
    name = newName;
}

std::string Material::getName() const
{
    return name;
}

void Material::setElementTagId(const int& value)
{
    elementTagId = value;
}

int Material::getElementTagId() const
{
    return elementTagId;
}


bool Material::readMaterialInputs()
{
    spdlog::debug("Material::readMaterialInputs");
    bool result {false};

    try
    {
        const std::string solverInputJsonFile = _directoryInputs + _materialInputFileName;
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

        if (materialJson.isArray())
        {
            spdlog::debug("Size of materialProperty array = {}", materialJson.size());

            for (int i = 0, n = materialJson.size(); i < n; i++)
            {
                readMaterialProperty(materialJson[i]);
            }
        }
        else 
        {
            // std::runtime_error("[ERROR] Solid::readMaterialInputs solverConfig.json file does't contain materialProperty array." );
            spdlog::info("materialProperty is not an array, single material defined");
            readMaterialProperty(materialJson);
        }

        result = true;
    }
    catch(const std::runtime_error& e)
    {
        std::cerr<<"[ERROR] Material::readMaterialInputs runtime exception = "<<e.what()<<"\n";
        spdlog::critical("Material::readMaterialInputs runtime exception = {}", e.what());
    }
    catch(const std::exception& e)
    {
        std::cerr<<"[ERROR] Material::readMaterialInputs std exception = " << e.what() << '\n';
        spdlog::critical("Material::readMaterialInputs std exception = {}", e.what());
    }
    catch (const char* errMsg)
    {
        std::cerr<<"Exception raised = "<<errMsg<<std::endl;
        spdlog::critical("Exception raised = {}", errMsg);
    }

    return result;
}


void Material::readMaterialProperty(const Json::Value &materialJson)
{
    // Temporary Variable used for reading inputs
    std::string tmp;

    if (materialJson.isMember("category") && materialJson["category"].isString())
    {
        tmp = materialJson["category"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        spdlog::info("Material type = {}", tmp);
        if (tmp == "SOLID") {
            Solid::readSolidMaterialInputs(materialJson, solid);
            spdlog::info("Read material solid={} property from json.", materialJson["name"].asString());
        }
        else if (tmp == "THERMAL") {
           Thermal::readThermalMaterialInputs(materialJson, thermal);
           spdlog::info("Read material thermal={} property from json.", materialJson["name"].asString());
        }
        else
        {
            spdlog::error("Unknown material category={} given in materialProperty of solverConfig", tmp);
            throw "Currently only SOLID and THERMAL material category are supported.";
        }
    }
    else
    {
        spdlog::error("Category field in material property={} , isCategory string={}", materialJson.isMember("category"), materialJson.isString());
        std::runtime_error("Category field in material missing or not string , check logs.");
    }
}
