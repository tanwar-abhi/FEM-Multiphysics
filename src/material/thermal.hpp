#pragma once

#include <Eigen/Core>

#include <jsoncpp/json/value.h>
#include <jsoncpp/json/json.h>

#include <string>
#include <fstream>
#include <iostream>
#include <optional>

#include "materialBase.hpp"


enum class ThermalMaterialType
{
    Insulator,
    Conductor,
    Thermal,
    Isotropic,
    Anisotropic,
    Unknown
};

// Data structure to hold material properties of domain for heat transfer problem.
class MaterialThermal
{
    private:
    ThermalMaterialType _materialType = ThermalMaterialType::Unknown;

    public:
    // Default Constructor and Destructors
    MaterialThermal();
    MaterialThermal(const std::string& inputsDirectoryPath);
    ~MaterialThermal();

    // Thickness of domain (element) to be defined for 2D cartesian coordinate problem only {units : m}
    // NOTE -> Never set thickness = 0, by default value of thickenss = 1
    double thickness = 1;

    // Mass Density {SI units = kg/m3}
    double massDensity = 0;

    // Specific heat (c) {SI units in J/Kg °C}
    double specificHeat = 0;

    /* Thermal conductivity for isotropic material {SI unit = W/(m °C)}
    k would be a tensor value for anisotropic material 
    In 2D, k = [Kxx, Kxy; Kyx, Kyy]; 
    In 3D k = [Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz]
    */
    double thermalConductivity2D[4] = {0, 0, 0, 0};
    double thermalConductivity3D[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    // Static function can be used without object, reads all material properties for Heat Transfer defined in inputs file
    static std::vector<std::optional<MaterialThermal>> readMaterialInputs();

    // void setInputDirectory(const std::string& inputJsonDirectory);

    // std::string getInputDirectory();

    Eigen::MatrixXd getMaterialThermalConductivity(const int& dimension);
};