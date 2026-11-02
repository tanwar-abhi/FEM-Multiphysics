
#include "thermal.hpp"


MaterialThermal::MaterialThermal()
{
    std::cout<<"[INFO] intialized thermal material constructor\n";
}


MaterialThermal::~MaterialThermal()
{

}


Eigen::MatrixXd MaterialThermal::getMaterialThermalConductivity(const int& dimension)
{
    Eigen::MatrixXd thermalConductivityMatrix(dimension, dimension);

    if (dimension == 2){
        int nonZero = 0;
        for (int i = 0; i < 4; i++){
            if (thermalConductivity2D[i] != 0){
                nonZero++;
            }
        }

        // Isothermal material (one value of thermal conductivity)
        if (nonZero == 1){
            thermalConductivityMatrix.row(0) << thermalConductivity2D[0], 0;
            thermalConductivityMatrix.row(1) << 0, thermalConductivity2D[0];
        }
        // Anisothermal material
        else{
            thermalConductivityMatrix.row(0) << thermalConductivity2D[0], thermalConductivity2D[1];
            thermalConductivityMatrix.row(1) << thermalConductivity2D[2], thermalConductivity2D[3];
        }
    }
    else if (dimension == 3){
        // Isothermal material (one value of thermal conductivity)
        thermalConductivityMatrix.row(0) << thermalConductivity2D[0] , 0 , 0;
        thermalConductivityMatrix.row(1) << 0 , thermalConductivity2D[0] , 0;
        thermalConductivityMatrix.row(2) << 0 , 0 , thermalConductivity2D[0];
    }

    return thermalConductivityMatrix;
}
