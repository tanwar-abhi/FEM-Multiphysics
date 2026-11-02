#include "materialBase.hpp"

template<typename T>
MaterialBase<T>::MaterialBase()
{
}

template<typename T>
MaterialBase<T>::~MaterialBase()
{
}

template<typename T>
void MaterialBase<T>::setName(const std::string& value)
{
    name = value;
}

template<typename T>
std::string MaterialBase<T>::getName() const
{
    return name;
}

template<typename T>
void MaterialBase<T>::setElementTagId(const int& tagID)
{
    elementTagId = tagID;
}

template<typename T>
int MaterialBase<T>::getElementTagId() const
{
    return elementTagId;
}

template<typename T>
std::string MaterialBase<T>::getInputsDirectory() const
{
    return _directoryInputs;
}

template<typename T>
std::string MaterialBase<T>::getMaterialInputFileName() const
{
    return _materialInputFileName;
}
