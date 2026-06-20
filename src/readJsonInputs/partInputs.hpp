#pragma once

#include <jsoncpp/json/json.h>
#include <jsoncpp/json/value.h>


#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <memory>
#include <vector>

class Part {
    public:
    Part();
    ~Part();
    Part(const std::string& m_name, const std::string& m_format, const std::string& m_meshFileName);

    std::string name;
    std::string format;
    std::string meshFileName;
};


class PartsInput {
    private:
    int numParts{0};

    const std::string fileName = "partsConfig.json";

    public:
    PartsInput();
    ~PartsInput();

    void setNumberOfParts(int);
    const int getNumberOfParts();

    // Deep Copy Constructor for PartsInp
    // PartsInput(const PartsInput &obj);

    // Part* part = NULL;
    std::vector<std::shared_ptr<Part>> parts;

    void readPartInputs(const std::string&);
};

