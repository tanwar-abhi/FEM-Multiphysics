#pragma once

#include <jsoncpp/json/value.h>
#include <jsoncpp/json/json.h>

#include <string>

class MaterialBase
{
    public:
    MaterialBase();
    ~MaterialBase();

    void setName(const std::string& name);
    std::string getName();

    void setElementTagId(const int& elementID);
    int getElementTagId();

    template <typename T>
    std::vector<T> readMaterialInputs();

    private:
    // Material Defination identifier name
    std::string name;

    // Element tag id to map material with equation and boundary
    int elementTagId;
};