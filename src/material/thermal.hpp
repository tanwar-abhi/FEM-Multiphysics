#pragma once

#include <jsoncpp/json/value.h>
#include <jsoncpp/json/json.h>

#include <string>

class Materials
{
    public:
    Materials();
    ~Materials();

    void setName(const std::string& name);
    std::string getName();

    void setElementTagId(const int& elementID);
    int getElementTagId();

    private:
    // Material Defination identifier name
    std::string name;

    // Element tag id to map material with equation and boundary
    int elementTagId;
};