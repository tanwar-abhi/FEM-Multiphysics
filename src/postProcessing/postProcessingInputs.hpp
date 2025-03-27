#include <fstream>
#include <json/value.h>
#include <json/json.h>
#include <iostream>
#include <string>

class PostProcess
{
    public:
    
    std::string type;

    int frequency;

    void readpostProcessingInputs();
};