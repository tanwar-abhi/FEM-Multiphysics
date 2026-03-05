#pragma once

#include <jsoncpp/json/json.h>
#include <jsoncpp/json/value.h>


#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

class Part {
    public:
    Part();
    ~Part();
    Part(const std::string m_name, const std::string m_format, const std::string m_meshFileName);

    std::string name;
    std::string format;
    std::string meshFileName;
};


class PartsInput {
    public:
    PartsInput();
    ~PartsInput();

    // Deep Copy Constructor for PartsInp
    PartsInput(const PartsInput &obj);

    int numParts = 0;
    Part* part = NULL;

    void readPartInputs();
};

