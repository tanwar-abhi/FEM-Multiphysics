
#include "partInputs.hpp"

#include <filesystem>

Part::Part()
{

}

Part::~Part()
{

}

Part::Part(const std::string& m_name, const std::string& m_format, const std::string& m_meshFileName)
{
    name = m_name;
    format = m_format;
    meshFileName = m_meshFileName;
}

PartsInput::PartsInput()
{

}

// Copy(Deep) constructor for PartsInput
PartsInput::PartsInput(const PartsInput &obj)
{
    numParts = obj.numParts;
    parts = obj.parts;
}

PartsInput::~PartsInput()
{

}

void PartsInput::setNumberOfParts(int value)
{
    numParts = value;
}

int PartsInput::getNumberOfParts()
{
    return numParts;
}

void PartsInput::readPartInputs(const std::string& filePath)
{
    std::cout<<"current filepath = " << filePath + "parts.json" <<std::endl;

    std::ifstream partsJsonFile(filePath + "parts.json");
    Json::Reader partsReader;
    Json::Value partsRoot;

    if(!partsReader.parse(partsJsonFile, partsRoot)){
        std::cerr << "Error while parsing parts"<<std::endl;
        exit(-404);
    }

    const Json::Value allParts = partsRoot["parts"];

    // Temporary string variable to be read into from PartInputs
    std::string tmp1;

    setNumberOfParts(allParts.size());

    for(int index = 0, n = getNumberOfParts(); index < n; index++)
    {
        std::shared_ptr<Part> part = std::make_shared<Part>();

        tmp1 = allParts[index]["format"].asString();
        std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::toupper);
        part->format = tmp1;

        tmp1 = allParts[index]["name"].asString();
        std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::toupper);
        part->name = tmp1;

        tmp1 = allParts[index]["meshFileName"].asString();
        part->meshFileName = tmp1;

        parts.emplace_back(part);
    }
}
