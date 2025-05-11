// Header file for read mesh functions


// NOTE:: This function only works on exported MeshFormat = 2.2 0 8
//        The mesh file must also contain node PhysicalNames, don't forget to define
//        the boundary PhysicalNames within GMSH.

#pragma once

#include "partInputs.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>

class readMesh{
    private:

    // Declaration of indigenous functions of readMesh, called by other functions within readmesh
    bool readElementTags(int &, int &, std::string);
    void pushNodeCoordString(std::string, std::vector<std::vector<double>> &, const int);
    void split(std::string Line, std::string &elemTagName, int &elemTagId);
    
    public:
    // Object of PartsInp class that is pases as input argument for readMesh constructor
    PartsInp partsInpObj;

    // Total number of elements
    // unsigned long int Nelements;

    // Total number of nodes in mesh
    unsigned long int NNodes;
    std::vector<unsigned long int> startingNodeIdMesh;

    // Node coordinates
    std::vector<std::vector<double>> Node_Coord;

    // Number of elemental tags
    int NElementalTags;

    // Physical field - elemental tags
    std::vector<std::string> elementalTags;
    std::vector<std::string> elementalTagsMeshName;
    std::vector<int> elementalTagsPartsId;
    std::vector<int> elementalTagIds;

    // Element Type number as per gmsh documentiation : https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    // BeType : Mesh bounday element type
    int elemType, BeType;

    // Store physical names tag of boundary nodes to identify during reading and interpolation.
    // order = element interpolation order; Nueman = NBtag; Dirichlet = DBtag;
    // int NBtag, DBtag;

    // Neuman boundary element connectivities, since GMSH considers boundary comprises of 1D elements
    // std::vector<std::vector<unsigned int>> NBConnect;

    // Dimension of each node
    int NODE_DIM;

    readMesh();
    ~readMesh();

    // Surfaces Boundary;

    // Constructor for readmesh object.
    readMesh(const PartsInp &partsInpObj, int dimension);

    // Read nodes coordinates from the msh file
    void readNodes();

    // Read Mesh size from msh file
    // void meshSize();

    // Return element connectivities and read Boundary nodes into respective set containers.
    void readElementTagsNew();
};


// Declaration of indigenous functions of readMesh (functional programing paradigm), called by other functions within readmesh
int countWords(std::string );

// Read the nodes on different boundaries defined in the mesh
// std::vector<std::vector<unsigned int>> ReadBoundaryNodes(std::string fileN, std::string boundaryName);

// Read the nodes on different boundaries defined in the mesh
std::vector<std::vector<unsigned int>> ReadBoundaryNodes(std::string fileN, std::string boundaryName);






// ################### Testing functions used only for debugging ##########################
// Print functions to check and debugg the read data.
// print 2D Vectors on terminal used for check and debugging
template<typename T>
void print2DVector(std::vector<std::vector<T>> toPrint)
{
    for (int i = 0, n = toPrint.size(); i < n; i++)
    {
        for (int j = 0, m = toPrint[i].size(); j < m; j++)
        {
            std::cout<<toPrint[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
}


// Print Boundary nodes from 1D vector
template<typename T>
void print1DVector(std::vector<T> myVector)
{
    for (int i = 0; i < myVector.size(); i++)
    {
        std::cout<<myVector[i]<<" ";
    }
    std::cout<<std::endl;
}
