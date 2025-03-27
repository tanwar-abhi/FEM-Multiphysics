#pragma once

#include <fstream>
#include <json/json.h>
#include <json/value.h>
#include <iostream>
#include <string>
#include <algorithm>

class Part{
    public:
    Part();
    ~Part();

    std::string name;
    std::string format;
    std::string meshFileName;
};


class PartsInp{
    public:
    PartsInp();
    ~PartsInp();

    // Deep Copy Constructor for PartsInp
    PartsInp(const PartsInp &obj);

    int numParts = 0;
    Part* part = NULL;

    void readPartInputs();
};

