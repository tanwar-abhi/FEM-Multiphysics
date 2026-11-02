#pragma once

#include <jsoncpp/json/value.h>
#include <jsoncpp/json/json.h>

#include <string>
#include <vector>

template <typename T>
class MaterialBase
{
    public:
    MaterialBase();
    virtual ~MaterialBase();
    MaterialBase(const std::string& inputsDirectory) : _directoryInputs{inputsDirectory}{}

    void setName(const std::string& name);
    std::string getName() const;

    void setElementTagId(const int& elementID);
    int getElementTagId() const;

    std::string getInputsDirectory() const;

    std::string getMaterialInputFileName() const;

    // virtual std::vector<T> readMaterialInputs() = 0;

    private:
    // Directory containing input json files full path
    std::string _directoryInputs;

    // Material Defination identifier name
    std::string name;

    // Element tag id to map material with equation and boundary
    int elementTagId;

    // file name that contains material inputs
    const std::string _materialInputFileName = "solverConfig.json";
};