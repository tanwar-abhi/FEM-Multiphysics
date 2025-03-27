/*
 * Element.cpp
 *
 *  Created on: 03-Oct-2021
 *      Author: pardhagv
 */
#include "Element.hpp"
#include <iostream>
#include <sstream>
#include <fstream>

Element::Element(){

}

Element::~Element(){
    // std::cout<<" # closure of element class obj #\n";
}

void Element::readElement(const readMesh &mesh, int tagId)
{
	int elemTagId = mesh.elementalTagIds[tagId];

	// This is useful for case when we have multiple seperate parts
	unsigned int startingNodeId = mesh.startingNodeIdMesh[mesh.elementalTagsPartsId[tagId]];

	// Open and read file.
	std::ifstream meshFile(mesh.elementalTagsMeshName[tagId]);

	// Read each line
	std::string Line;

	// Element connectivities reading flag
	bool eReadFlag = false;

	// flag for checking if node already exists
	int flag;

	//
	// int nPosition;

	int tagColId;

	this->numElems = 0;

	while (getline(meshFile, Line))
	{
		// Removing any empty space or carraige reurn operators from line
		if (!Line.empty() && Line[Line.size() - 1] == '\r')
			Line.erase(Line.size() - 1);

		if (Line == "$Elements")
		{
			eReadFlag = true;
			continue;
		}

		if (Line == "$EndElements")
		{
			meshFile.close();

			// Determining the number of nodes and elminating duplicate nodes
			std::sort(this->Nodes.begin(), this->Nodes.end());
			auto last = std::unique(this->Nodes.begin(), this->Nodes.end());
			this->Nodes.erase(last, this->Nodes.end());
			this->numNodes = this->Nodes.size();

			// Debugging part
			// std::cout<<"Number of nodes "<<this->numNodes<<std::endl;
			// std::cout<<"Number of elements "<<this->numElems<<std::endl;

			// ###################################
			// lines for debugging
			// std::cout<<"Connectivity = \n";
			// print2DVector(this->ElemConnectivity);
			// std::cout<<"Nodes = \n";
			// print1DVector(this->Nodes);
      		// std::cout<<std::endl;
			// ###################################

			// std::cout<<"\n # eType = "<<this->elemType<<" "<<this->numElems<<" "<<this->numNodes <<std::endl;

			// ###################################################

			return;
		}

		// Count the number of columns in the line
		// If the number of columns is 8
		// Three noded triangle and second order line elements
		if (eReadFlag)
		{
			// One dimensional line element with 2 nodes
			if (countWords(Line) ==  7)
			{
				std::stringstream SS(Line);

				// Element Type (etype), total no. of tags (tagQty), Element Tag (etag) as per gmsh file format
				int _, etype, tagQty, etag;
				unsigned int node1, node2;
				std::vector<unsigned int> lclElemCnn;
				SS  >> _ >> etype >> tagQty >> etag;
				// std::vector<unsigned int> nbT;

				// Read Redundant/unnecessary tags
				for (int i = 1; i < tagQty; i++) {SS >> _;}

				if (etag == elemTagId)
				{
					// 2 node line
					if (etype == 1) 
					{
						SS >> node1 >> node2;
						// elemental nodes
						lclElemCnn.push_back(node1 + startingNodeId);
						lclElemCnn.push_back(node2 + startingNodeId);
						// Nodes
						this->Nodes.push_back(node1 + startingNodeId);
						this->Nodes.push_back(node2 + startingNodeId);
						
						this->numNodesElem = 2;
						this->basisFn1D.getShapeFn(etype);
					}
					this->ElemConnectivity.push_back(lclElemCnn);
					this->numElems++;
					this->elemType = etype;
				}
			}

			if (countWords(Line) ==  8)
			{
				std::stringstream SS(Line);
				int _, etype, tagQty, etag;
				unsigned int node1, node2, node3;
				std::vector<unsigned int> lclElemCnn;
				SS  >> _ >> etype >> tagQty >> etag ;

				// Read Redundant/unnecessary tags
				for (int i = 1; i < tagQty; i++) {SS >> _;}

				// 3 noded linear triange or 3 noded line element
				if (etag == elemTagId)
				{

					// Three noded triangle 2D element
					if (etype == 2) 
					{
						SS >> node1 >> node2 >> node3;
						lclElemCnn.push_back(node1+startingNodeId);
						lclElemCnn.push_back(node2+startingNodeId);
						lclElemCnn.push_back(node3+startingNodeId);
						// Nodes
						this->Nodes.push_back(node1 + startingNodeId);
						this->Nodes.push_back(node2 + startingNodeId);
						this->Nodes.push_back(node3 + startingNodeId);

						this->numNodesElem = 3;
						this->basisFn2D.getShapeFn(etype, 3);
					}
					// Three noded second order line
					else if (etype == 8)
					{
						SS >> node1 >> node2 >> node3;
						lclElemCnn.push_back(node1 + startingNodeId);
						lclElemCnn.push_back(node2 + startingNodeId);
						lclElemCnn.push_back(node3 + startingNodeId);
						// Nodes
						this->Nodes.push_back(node1 + startingNodeId);
						this->Nodes.push_back(node2 + startingNodeId);
						this->Nodes.push_back(node3 + startingNodeId);

						this->numNodesElem = 3;
						this->basisFn1D.getShapeFn(etype, 3);
					}
					this->elemType = etype;
					this->ElemConnectivity.push_back(lclElemCnn);
					this->numElems++;
				}
			}
			// 4 Noded Quad element or 4 Noded tetrahedron
			else if (countWords(Line) ==  9) 
			{
				std::stringstream SS(Line);
				int _, etype, tagQty, etag;
				unsigned int node1, node2, node3, node4;
				std::vector<unsigned int> lclElemCnn;
				SS  >> _ >> etype >> tagQty >> etag ;

				// Read Redundant tags
				for (int i = 1; i < tagQty; i++) {SS >> _;}

				if (etag == elemTagId)
				{
					// 4 noded, 2D Quad
					if (etype == 3) 
					{
						SS >> node1 >> node2 >> node3 >> node4;
						lclElemCnn.push_back(node1 + startingNodeId);
						lclElemCnn.push_back(node2 + startingNodeId);
						lclElemCnn.push_back(node3 + startingNodeId);
						lclElemCnn.push_back(node4 + startingNodeId);
						// Nodes
						this->Nodes.push_back(node1 + startingNodeId);
						this->Nodes.push_back(node2 + startingNodeId);
						this->Nodes.push_back(node3 + startingNodeId);
						this->Nodes.push_back(node4 + startingNodeId);

						this->numNodesElem = 4;
						this->basisFn2D.getShapeFn(etype);
					}

					// 4 node tetrahedron element
					else if (etype == 4)
					{
						// ####### TODO ###########
					}
					this->elemType = etype;
					this->ElemConnectivity.push_back(lclElemCnn);
					this->numElems++;
				}
			}
			// 6 Noded Triangular element
			else if (countWords(Line) ==  11)
			{
				std::stringstream SS(Line);
				int _, etype, tagQty, etag;
				unsigned int node1, node2, node3, node4, node5, node6;
				std::vector<unsigned int> lclElemCnn;
				SS >> _ >> etype >> tagQty >> etag;

				// Read Redundant/unnecessary tags
				for (int i = 1; i < tagQty; i++) {SS >> _;}

				if (etag == elemTagId)
				{
					// 6 Noded Triangular element
					if (etype == 9)
					{
						SS >> node1 >> node2 >> node3 >> node4 >> node5 >> node6;
						lclElemCnn.push_back(node1 + startingNodeId);
						lclElemCnn.push_back(node2 + startingNodeId);
						lclElemCnn.push_back(node3 + startingNodeId);
						lclElemCnn.push_back(node4 + startingNodeId);
						lclElemCnn.push_back(node5 + startingNodeId);
						lclElemCnn.push_back(node6 + startingNodeId);
						// Nodes
						this->Nodes.push_back(node1 + startingNodeId);
						this->Nodes.push_back(node2 + startingNodeId);
						this->Nodes.push_back(node3 + startingNodeId);
						this->Nodes.push_back(node4 + startingNodeId);
						this->Nodes.push_back(node5 + startingNodeId);
						this->Nodes.push_back(node6 + startingNodeId);

						this->numNodesElem = 6;
						this->basisFn2D.getShapeFn(etype);
					}
					this-> elemType = etype;
					this->ElemConnectivity.push_back(lclElemCnn);
					this->numElems++;
				}
			}
			// 8 noded hexahedral element
			else if (countWords(Line) == 13)
			{
				std::stringstream ss(Line);
				int _, etype, tagQty, etag;
				// unsigned int node1, node2, node3, node4, node5, node6, node7, node8;
				unsigned int nodes;
				std::vector<unsigned int> lclElemCnn;
				ss >> _ >> etype >> tagQty >> etag;

				// Read Redundant/unnecessary tags
				for (int i = 1; i < tagQty; i++) {ss >> _ ;}

				if (etag == elemTagId)
				{
					// 8 noded hexahedral element {3D element}
					if (etype == 5)
					{
						this->numNodesElem = 8;
						this->basisFn3D.getShapeFn(etype, 8);
						for (int i = 0; i < numNodesElem; i++)
						{
							ss >> nodes;
							lclElemCnn.push_back(nodes + startingNodeId);
							this->Nodes.push_back(nodes + startingNodeId);
						}
					}
					this->elemType = etype;
					this->ElemConnectivity.push_back(lclElemCnn);
					this->numElems++;
				}
			}
			// 9 Noded Quad element
			else if (countWords(Line) ==  14)
			{
				std::stringstream SS(Line);
				int _, etype, tagQty, etag;
				unsigned int node1, node2, node3, node4, node5, node6, node7, node8, node9;
				std::vector<unsigned int> lclElemCnn;
				SS >> _ >> etype >> tagQty >> etag;
				// std::vector<unsigned int> nbT;

				// Read Redundant/unnecessary tags
				for (int i = 1; i < tagQty; i++) {SS >> _;}

				// 9 Noded Quad element {second order element}
				if (etag == elemTagId)
				{
					if (etype == 10)
					{
						SS >> node1 >> node2 >> node3 >> node4 >> node5 >> node6 >> node7 >> node8 >> node9;
						lclElemCnn.push_back(node1 + startingNodeId);
						lclElemCnn.push_back(node2 + startingNodeId);
						lclElemCnn.push_back(node3 + startingNodeId);
						lclElemCnn.push_back(node4 + startingNodeId);
						lclElemCnn.push_back(node5 + startingNodeId);
						lclElemCnn.push_back(node6 + startingNodeId);
						lclElemCnn.push_back(node7 + startingNodeId);
						lclElemCnn.push_back(node8 + startingNodeId);
						lclElemCnn.push_back(node9 + startingNodeId);
						// Nodes
						this->Nodes.push_back(node1 + startingNodeId);
						this->Nodes.push_back(node2 + startingNodeId);
						this->Nodes.push_back(node3 + startingNodeId);
						this->Nodes.push_back(node4 + startingNodeId);
						this->Nodes.push_back(node5 + startingNodeId);
						this->Nodes.push_back(node6 + startingNodeId);
						this->Nodes.push_back(node7 + startingNodeId);
						this->Nodes.push_back(node8 + startingNodeId);
						this->Nodes.push_back(node9 + startingNodeId);

						this->elemType = etype;
						this->numNodesElem = 9;
						this->ElemConnectivity.push_back(lclElemCnn);
						this->numElems++;

						this->basisFn2D.getShapeFn(etype);
					}
				}
			}
		}
	}
}

