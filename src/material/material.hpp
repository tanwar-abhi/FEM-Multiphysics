
#pragma once

#include "solid.hpp"
#include "thermal.hpp"

#include <optional>


// class Material : public Solid, public MaterialThermal 
class Material
{
	public:
		Material();
		~Material();

		Material(const std::string& inputDirectory) : _directoryInputs{inputDirectory} {};

		std::vector<Solid> solid;
		std::vector<Thermal> thermal;

		bool hasSolid() const;
		bool hasThermal() const;

		void setName(const std::string& name);
		std::string getName() const;

		void setElementTagId(const int& elementID);
		int getElementTagId() const;

		std::string getInputsDirectory() const;
		void setInputsDirectory(const std::string& inputsDir);

		std::string getMaterialInputFileName() const;

		bool readMaterialInputs();

	private:
		// Directory containing input json files full path
		std::string _directoryInputs;

		// Material Defination identifier name
		std::string name;

		// Element tag id to map material with equation and boundary
		int elementTagId;

		// file name that contains material inputs
		const std::string _materialInputFileName = "solverConfig.json";

		void readMaterialProperty(const Json::Value &materialJson);

};
