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

std::string MaterialBase::getName() const
{
    return name;
}

void MaterialBase::setElementTagId(const int& tagID)
{
    elementTagId = tagID;
}

int MaterialBase::getElementTagId() const
{
    return elementTagId;
}

void setInputsDirectory(const std::string& inputsDirectroyWithFullPath)
{
    _directoryInputs = inputsDirectoryWithFullPath;
}

std::string getInputDirectory() const
{
    return _directoryInputs;
}
