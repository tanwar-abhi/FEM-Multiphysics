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
    return solid[0].has_value();
}

bool Material::hasThermal() const
{
	return thermal[0].has_value();
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
    bool result {false};
    try
    {
        Solid::readSolidMaterialInputs(getInputsDirectory() + getMaterialInputFileName(), solid);
        spdlog::debug("Material::readMaterialInputs readSolidMaterialInputs read material config success");

        // thermal = MaterialThermal::readMaterialInputs();
        // spdlog::debug("Material::readMaterialInputs readMaterialInputs read material config success");
        result = true;
    }
    catch(const std::runtime_error& e)
    {
        std::cerr<<"[ERROR] Material::readMaterialInputs runtime exception = "<<e.what()<<"\n";
    }
    catch(const std::exception& e)
    {
        std::cerr<<"[ERROR] Material::readMaterialInputs std exception = " << e.what() << '\n';
    }
    
    return result;
}




