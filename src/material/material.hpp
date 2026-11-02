
#pragma once

#include <solid.hpp>
#include <thermal.hpp>

#include <optional>


class Material : public Solid, public MaterialThermal 
{
	public:
		Material();
		~Material();

		std::optional<Solid> solid;
		std::optional<MaterialThermal> thermal;

		bool hasSolid() const;
		bool hasThermal() const;

		void setName(const std::string& name);
		std::string getName() const;
		
		void setElementTagId(const int& elementID);
		int getElementTagId() const;
		
		std::string getInputsDirectory() const;
		
		std::string getMaterialInputFileName() const;
		
		std::vector<Material> readMaterialInputs();

	private:
		// Directory containing input json files full path
		std::string _directoryInputs;

		// Material Defination identifier name
		std::string name;

		// Element tag id to map material with equation and boundary
		int elementTagId;

		// file name that contains material inputs
		const std::string _materialInputFileName = "solverConfig.json";

};
