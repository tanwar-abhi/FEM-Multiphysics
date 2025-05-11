
#include "partInputs.hpp"

Part::Part()
{

}

Part::~Part()
{

}

Part::Part(const std::string m_name, const std::string m_format, const std::string m_meshFileName)
{
    name = m_name;
    format = m_format;
    meshFileName = m_meshFileName;
}

PartsInp::PartsInp()
{

}

// Copy(Deep) constructor for partsInp, to be called whenever an object is called as argument to a function.
PartsInp::PartsInp(const PartsInp &obj)
{
    numParts = obj.numParts;

    // Allocate new memory from Heap for copy
    part = new Part[numParts];

    for (int i = 0; i < numParts; i++)
    {
        part[i] = obj.part[i];
    }
}


PartsInp::~PartsInp()
{

    // Release the block of memory
    // std::cout<<" # Part Input Destructor #"<<std::endl;
    delete[] part;
}

void PartsInp::readPartInputs()
{
    std::ifstream part_text ("parts.json");
    Json::Reader part_reader;
    Json::Value part_root;

    //test to check if parsing is happening properly
    bool parsingSuccessful = part_reader.parse(part_text, part_root);
    if(!parsingSuccessful){
        std::cerr << "Error while parsing parts"<<std::endl;
        exit(-404);
    }

    const Json::Value inp_parts = part_root["parts"];

    // Temporary string variable to be read into from PartInputs
    std::string tmp1;

    this->numParts = inp_parts.size();
    Part *part = new Part[inp_parts.size()];
    this->part = part;

    for(int index = 0; index < this->numParts; index++){
        tmp1 = inp_parts[index]["format"].asString();
        std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::toupper);

        this->part[index].format = tmp1;

        tmp1 = inp_parts[index]["name"].asString();
        std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::toupper);

        this->part[index].name = tmp1;

        tmp1 = inp_parts[index]["meshFileName"].asString();

        this->part[index].meshFileName = tmp1;
    }
}
