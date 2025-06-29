#include <fstream>
#include <jsoncpp/json/value.h>
#include <jsoncpp/json/json.h>
#include <iostream>
#include <string>

class PostProcess
{
    public:
    
    std::string type;

    int frequency;

    void readpostProcessingInputs();
};