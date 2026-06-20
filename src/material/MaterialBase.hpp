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
    ~MaterialBase();

    void setName(const std::string& name);
    std::string getName() const;

    void setElementTagId(const int& elementID);
    int getElementTagId() const;

    void setInputsDirectory(const std::string&  inputDirectory);
    std::string getInputsDirectory() const;

    virtual std::vector<T> readMaterialInputs() = 0;

    private:
    // Directory containing input json files full path
    static std::string _directoryInputs;

    // Material Defination identifier name
    std::string name;

    // Element tag id to map material with equation and boundary
    int elementTagId;

    // Directory containing material inputs file with whole path
    std::string _directoryInputs;
};