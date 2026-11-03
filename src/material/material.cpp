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

std::string Material::getInputsDirectory() const
{
    return _directoryInputs;
}

bool Material::readMaterialInputs()
{
    bool result {false};
    try
    {
        solid = Solid::readMaterialInputs();
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




