#include "postProcessingInputs.hpp"

PostProcess postProcess;

void readpostProcessInputs(){
    std::ifstream text("postProcessingInputs");
    Json::Value postProcess_root;
    Json::Reader postProcess_reader;
    bool parsingSuccessful = postProcess_reader.parse(text, postProcess_root);
    if(!parsingSuccessful)
    {
	std::cout << "Error while parsing" << std::endl;        
    }

    const Json::Value inp_postProcessing = postProcess_root["nodalOutputs"];

    postProcess.type = inp_postProcessing[0]["type"].asString();
    postProcess.frequency = inp_postProcessing[0]["frequency"].asInt();

}