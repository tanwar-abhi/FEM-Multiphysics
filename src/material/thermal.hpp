#pragma once

#include <jsoncpp/json/value.h>
#include <jsoncpp/json/json.h>

#include <string>

// Data structure to hold material properties of domain for heat transfer problem.
class MaterialThermal : public MaterialBase
{
    public:
    // Default Constructor and Destructors
    MaterialThermal();
    ~MaterialThermal();

    // Thickness of domain (element) to be defined for 2D cartesian coordinate problem only {units : m}
    // NOTE -> Never set thickness = 0, by default value of thickenss = 1
    double thickness = 1;


    // Mass Density {SI units = kg/m3}
    double RHO = 0;

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
    static std::vector<MaterialThermal> readMaterialInputs();

    void setInputDirectory(const std::string& inputJsonDirectory);

    std::string getInputDirectory();
};