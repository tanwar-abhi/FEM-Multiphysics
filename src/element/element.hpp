/*
 * element.hpp
 *
 *  Created on: 03-Oct-2021
 *      Author: pardhagv
*/

#ifndef SRC_ELEMENT_HPP_
#define SRC_ELEMENT_HPP_

#include "readMesh.hpp"
#include "basisFunction.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <set>

class Element{
    private:

    public:
	// constructor and destructors
	Element();
	~Element();
	void readElement(const readMesh &, int );
 
	// Total number of elements
	int numElems;

	// Number of nodes for each element
	int numNodesElem;

	// Total number of nodes
	int numNodes;

	// Element type and order
	int elemType;
	int elemOrder;

	// Connectivity of elements
	std::vector<std::vector<unsigned int>> ElemConnectivity;

	// All (unique) nodes on the selected boundary in sorted order
	std::vector<unsigned int> Nodes;

	// Shape function objects for respective element
	ShapeFn2D basisFn2D;
	ShapeFn1D basisFn1D;
	ShapeFn3D basisFn3D;
};

void checkNodes(int* , std::vector<int> , int , int* );
#endif /* SRC_ELEMENT_HPP_ */
