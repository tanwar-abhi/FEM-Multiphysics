
#include "neumanBC.hpp"


void applyThermalNBC(std::vector<Tr> &tripletStiffness, Eigen::VectorXd &f_Global, const BoundaryConditions &boundary, 
                        const Element meshElements[], const readMesh &mesh, const SolverInp &solverInpObj, const MaterialThermal &material)
{

    // Apply neumann Boundary conditions
    for (int index = 0; index < boundary.nNBC; index++)
    {

        // If Surface Heat flux defined (at conduction boundary) then calculate contribution to heat vector
        if (boundary.neumann[index].variable == "HEATFLUX")
        {
            // Read nodes on the boundary where surface heating is defined
            int elementTag = boundary.neumann[index].elemTagId;
            // Heat load vector contribution from Surface heating over surface area (S2)
            SurfaceConduction(f_Global, mesh.Node_Coord, meshElements[elementTag], solverInpObj, boundary.neumann[index], material);
        }

        // Adding convection boundary contributions
        if (boundary.neumann[index].variable == "CONVECTIVEHEATTRANSFER")
        {
            if (boundary.neumann[index].H != 0)
            {
                // Get convection boundary nodes tag from ELEMENT
                int elementTag = boundary.neumann[index].elemTagId;

                // Contribution to the stiffness matrix due to convection boundary
                StiffnessConvection(tripletStiffness, mesh.Node_Coord, meshElements[elementTag], solverInpObj, boundary.neumann[index], material);

                // Surface convection term contribution to thermal load vector
                SurfaceConvection(f_Global, mesh.Node_Coord, meshElements[elementTag], solverInpObj, boundary.neumann[index], material);
            }
        }

        /* Line Heat Source/Sink boundary
        if (boundary.neuman[index].variable.find("SOURCE") != std::string::npos || boundary.neuman[index].variable.find("SINK") != std::string::npos)
        {
            if (boundary.neuman[index].Q != 0)
            {

                std::vector<std::vector<unsigned int>> heatSourceNodes = meshElements[boundary.neuman[index].elemTagId].ElemConnectivity;
                // std::vector<std::vector<double>> NODE_COORD = meshElements[boundary.neuman[index].elemTagId].Nodes;
                // heatSourceNodes = ReadBoundaryNodes(mesh.fileName, "source");

                LineHeatSource(f_Global, mesh.Node_Coord, mesh, heatSourceNodes, solverInpObj, material_Thermal[0], boundary.neuman[index]);
            }
        }*/
    }
}





void applyStructureNBC(Eigen::VectorXd &f, const readMesh &mesh, const Element meshElements[], const SolverInp &solverInpObj, const LinearElasticMaterial &material_LE, const BoundaryConditions &boundary, int eqnIndex)
{
    for (int i = 0; i<boundary.nNBC; i++)
    {
        if (boundary.neumann[i].variable == "TRACTION" || boundary.neumann[i].variable == "PRESSURE" || boundary.neumann[i].variable == "FORCE")
        {
            int tractionElemTag = boundary.neumann[i].elemTagId;
            imposeTraction(f, meshElements[tractionElemTag], mesh, solverInpObj, boundary.neumann[i], material_LE, eqnIndex);
        }
    }
}

