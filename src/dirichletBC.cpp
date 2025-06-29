
#include "dirichletBC.hpp"


// Search the DB induced node number data to obtain rows and columns number to be eliminated from stiffness matrix
bool SearchIn(std::map<unsigned int, double> &nodeValuePair, unsigned int row, unsigned int col)
{
	std::map<unsigned int, double>::iterator posRow = nodeValuePair.find(row);
	std::map<unsigned int, double>::iterator posCol = nodeValuePair.find(col);

	if (posRow != nodeValuePair.end() || posCol != nodeValuePair.end())
	{
		// Row or column index Found in vector
		return true;
	}

	// Not found in vector of dirichlet nodes
	return false;
}




// Function to get the new row and column numbers after elimination
int Counter(std::map<unsigned int, double> &nodePair, unsigned int index)
{
	int count = 0;

	for (std::map<unsigned int, double>::iterator itr = nodePair.begin(); itr != nodePair.end(); itr++)
	{
		if (itr->first < index)
			count++;
		else{
			break;
		}
	}
	return count;
}




// Get the respective node number as per DOF and values for vector field variables i.e. displacement or velocity
void getNodesPerDof(std::map<unsigned int, double> &nodeValuePair, const DirichletBC &boundary, unsigned int dbNode, double dirichletBvalues, int dof, int &j)
{
	if (boundary.variable == "DISPX" || boundary.variable == "VELX")
	{
		if (j == 0)
		{
			nodeValuePair.insert(std::pair<unsigned int, double>((dof * (dbNode - 1) + j), dirichletBvalues));
			j += dof;
		}
	}
	else if (boundary.variable == "DISPY" || boundary.variable == "VELY")
	{
		if (j == 1)
		{
			nodeValuePair.insert(std::pair<unsigned int, double>((dof * (dbNode - 1) + j), dirichletBvalues));
			j += dof;
		}
	}
	else if (boundary.variable == "DISPZ" || boundary.variable == "VELZ")
	{
		if (j == 2)
		{
			nodeValuePair.insert(std::pair<unsigned int, double>((dof * (dbNode - 1) + j), dirichletBvalues));
			j += dof;
		}
	}
}



/*/ ########### STABLE OLD method ##################################################
// Eliminate Dirichlet nodes from stiffness matrix (triplet sparse)
void EliminateDirichletNodes(std::vector<Eigen::Triplet<double>> &trList, Eigen::VectorXd &fglobal, std::map<unsigned int, double> &nodeValuePair)
{

	std::map<unsigned int, double>::iterator itr;
	for (itr = nodeValuePair.begin(); itr != nodeValuePair.end(); itr++)
	{
		// row or column has non zero values in stiffness matrix for given DB nodes
		for (int index = 0; index < trList.size(); index++)
		{
			if (trList[index].col() == itr->first)
			{
				// Move stiffness coefficient*node value to RHS vector
				fglobal(trList[index].row()) -= trList[index].value() * itr->second;

				// Remove column element that's been moved to RHS force vector
				trList.erase(trList.begin() + index);

				// Adjusting index and size of triplet after elimination
				--index;
			}

			// Reset negative index counter to check rows
			if (index < 0)
			{
				index = 0;
			}

			if (trList[index].row() == itr->first)
			{
				// Remove row element to reduce stiffness matrix
				trList.erase(trList.begin() + index);

				// Adjusting index and total size of triplet after elimination
				--index;
			}
		}
	}

	// Reset node row and column numbers after elimination of dirichlet nodes.
	for (int i = 0, n = trList.size(); i < n; i++)
	{
		// Index of row and column
		unsigned int iR = 0, iC = 0;
		
		if (!SearchIn(nodeValuePair, trList[i].row(), trList[i].col())){
			iR = Counter(nodeValuePair, trList[i].row());
			iC = Counter(nodeValuePair, trList[i].col());
			trList[i] = Eigen::Triplet<double>(trList[i].row() - iR, trList[i].col() - iC, trList[i].value());
		}
	}
}

// #########################################################################################
*/




// ######## Experimental new method to try to reduce time complexity #########
// Eliminate Dirichlet nodes from stiffness matrix (triplet sparse)
void EliminateDirichletNodes(std::vector<Eigen::Triplet<double>> &trList, Eigen::VectorXd &fglobal, std::map<unsigned int, double> &nodeValuePair)
{

	// std::vector<Eigen::Triplet<double>> tripletReduced;


	std::map<unsigned int, double>::iterator itr;
	for (itr = nodeValuePair.begin(); itr != nodeValuePair.end(); itr++)
	{
		// row or column has non zero values in stiffness matrix for given DB nodes
		for (int index = 0; index < trList.size(); index++)
		{
			if (trList[index].col() == itr->first)
			{
				// Move stiffness coefficient*node value to RHS vector
				fglobal(trList[index].row()) -= trList[index].value() * itr->second;

				// Remove column element that's been moved to RHS force vector
				trList.erase(trList.begin() + index);

				// Adjusting index and size of triplet after elimination
				--index;
			}

			// Reset negative index counter to check rows
			if (index < 0)
			{
				index = 0;
			}

			if (trList[index].row() == itr->first)
			{
				// Remove row element to reduce stiffness matrix
				trList.erase(trList.begin() + index);

				// Adjusting index and total size of triplet after elimination
				--index;
			}
		}
	}

	// Reset node row and column numbers after elimination of dirichlet nodes.
	for (int i = 0, n = trList.size(); i < n; i++)
	{
		// Index of row and column
		unsigned int iR = 0, iC = 0;
		
		if (!SearchIn(nodeValuePair, trList[i].row(), trList[i].col())){
			iR = Counter(nodeValuePair, trList[i].row());
			iC = Counter(nodeValuePair, trList[i].col());
			trList[i] = Eigen::Triplet<double>(trList[i].row() - iR, trList[i].col() - iC, trList[i].value());
		}
	}
}




// Elimnate the RHS vector (i.e. force vector) to impose Dirichlet boundary conditions
void eliminateRHS(Eigen::VectorXd &fglobal, std::map<unsigned int, double> &nodeValuePair)
{

	// Vector of domain node values of fglobal vector containg non DB nodes
	std::vector<double> fValues;

	for (unsigned int i = 0; i < fglobal.size(); ++i)
	{
        std::map<unsigned int, double>::iterator itr = nodeValuePair.find(i);

        if (itr == nodeValuePair.end())
        {
            fValues.push_back(fglobal(i));
        }
	}
	
    // Since resizing a vector is computationally expensive in Eigen thus created new vector
	// Initialise reduced (f_reduced) from fvalues and set it equal to fglobal to obtain eliminated force vector
	Eigen::Map<Eigen::VectorXd> f_reduced(fValues.data(), fValues.size());
	fglobal = f_reduced;
}




// Return solution after imposing Dirichlet boundary values at subsequesnt nodes
Eigen::MatrixXd ImposeDBC(const Eigen::VectorXd &u_red, std::map<unsigned int, double> &nodeValuePair, unsigned long numNodesDomain, int dof)
{

	// Solution vector
	Eigen::VectorXd result(numNodesDomain*dof);
    result.setZero();

	// Index counter for dirichlet nodes and result vector
	int count = 0;

    // Create vector of reduced value
    for (int i = 0; i < numNodesDomain*dof; i++)
    {
		// Iteration not in dirichlet nodes global value
		if (nodeValuePair.find(i) == nodeValuePair.end()){
			result(i) = u_red(count);
			count++;
		}
    }

    // Insert values using node value pair hash map
    std::map<unsigned int, double>::iterator itr;
    for (itr = nodeValuePair.begin(); itr != nodeValuePair.end(); itr++)
    {
        result(itr->first) = itr->second;
    }

	Eigen::MatrixXd solution(numNodesDomain, dof);

	if (dof == 1)
		return result;
	else
	{
		solution.setZero();

		// Index counter
		int ct = 0;

		// Map solution into size(NumNodes, dof), i.e. each row is number of nodes and column correspond to respective unknown (U) as per DOF
		// For example -> vector Problem Solution = [Ux, Uy, Uz] ; scalar Problem Solution = [U]
		for (int i = 0; i < numNodesDomain; i++){
			for (int j = 0; j < dof; j++){
				solution(i,j) = result(ct++);
			}
		}
	}

	return solution;
}




// void applyDirichletBC(std::vector<Eigen::Triplet<double>> &tripletGlobal, Eigen::VectorXd &fGlobal, const Element meshElement[], const BoundaryConditions &boundary, const SolverInp &solverInp)
Eigen::SparseMatrix<double> applyDirichletBC(std::vector<Eigen::Triplet<double>> &tripletGlobal, Eigen::VectorXd &fGlobal, const Element meshElement[], const BoundaryConditions &boundary, const SolverInp &solverInp, unsigned long numNodesDomain)
{
	// Total number of nodes(dirichlet) to be removed from global stiffness and force vector
	int numReducedNodes = 0;

	if (boundary.nDBC != 0)
	{
		// std::vector<unsigned int> dof(solverInp.nEquations);
		// ########### Considering single equation ########### to be changed in future
		unsigned int dof = solverInp.equations[0].DOF;

		// Vector of all dirichlet nodes and value at each dirichlet boundary
		std::vector<std::vector<unsigned int>> dirichletBNodes(boundary.nDBC);
		std::vector<double> dirichletBvalues(boundary.nDBC);

		for (int i = 0; i < boundary.nDBC; i++)
		{
			int dirichletTag = boundary.dirichlet[i].elemTagId;
			dirichletBNodes[i] = meshElement[dirichletTag].Nodes;
			dirichletBvalues[i] = boundary.dirichlet[i].value;
		}

		// Hash map (Key:Value pairs) of global dirichlet node number and value at those nodes
		std::map<unsigned int, double> nodeValuePair;


		for (int ie = 0, ieN = dirichletBNodes.size(); ie < ieN; ie++)
		{
			for (unsigned int dbNode : dirichletBNodes[ie])
			{
				for (int j = 0; j < dof; j++)
				{
					if (boundary.dirichlet[ie].variable.find("DISP") != std::string::npos || boundary.dirichlet[ie].variable.find("VEL") != std::string::npos){
						// In case of vector field problems, get the node and value to be prescribed at the node as per DOF
						getNodesPerDof(nodeValuePair, boundary.dirichlet[ie], dbNode, dirichletBvalues[ie], dof, j);
					}
					else{
						// Scalar field problem, i.e. single unknown per node (DOF = 1)
						nodeValuePair.insert(std::pair<unsigned int, double>((dof * (dbNode - 1) + j), dirichletBvalues[ie]));
					}
				}
			}
		}

		// EliminateDirichletNodes(tripletGlobal, fGlobal, e2dof, nodeValuePair);
		EliminateDirichletNodes(tripletGlobal, fGlobal, nodeValuePair);

		eliminateRHS(fGlobal, nodeValuePair);

		// Total number of nodes to be removed from stiffness matrix (dirichlet nodes as per DOF)
		numReducedNodes = nodeValuePair.size();
	}

	// Final size of stiffness matrix after elimination of dirichlet nodes
	unsigned long reducedSize = (numNodesDomain*solverInp.equations[0].DOF)-(numReducedNodes);

	// Reduced stiffness matrix after elimination
	Eigen::SparseMatrix<double> Kreduced(reducedSize, reducedSize);
	Kreduced.setFromTriplets(tripletGlobal.begin(), tripletGlobal.end());

	return Kreduced;
}






// Impose dirichlet boundary values at solution to satisfy dirichlet conditions
Eigen::MatrixXd imposeDirichletValues(const Eigen::VectorXd &solutionRed, const Element meshElement[], const BoundaryConditions &boundary, unsigned long numNodesDomain, int dof)
{
	// Vector of all dirichlet nodes
	std::vector<std::vector<unsigned int>> dirichletBNodes(boundary.nDBC);

	// Vecotor of all values at the subsequesnt boundary
	std::vector<double> dirichletBValues(boundary.nDBC);

	for (int i = 0; i <boundary.nDBC; i++ ){
		int dirichletTag = boundary.dirichlet[i].elemTagId;

		dirichletBNodes[i] = meshElement[dirichletTag].Nodes;
		dirichletBValues[i] = boundary.dirichlet[i].value;
	}


	//Hash map (Key:Value pairs) of global node number and vlaue at that node
	std::map<unsigned int, double> nodeValuePair;

	for (int i = 0, n = dirichletBNodes.size(); i<n; i++)
	{
		for (double dbNode : dirichletBNodes[i])
		{
			for (int j = 0; j < dof; j++){

				if (boundary.dirichlet[i].variable.find("DISP") != std::string::npos || boundary.dirichlet[i].variable.find("VEL") != std::string::npos){
					// In case of vector field problems, get the node and value to be prescribed at the node as per DOF
					getNodesPerDof(nodeValuePair, boundary.dirichlet[i], dbNode, dirichletBValues[i], dof, j);
				}
				else{
					// Scalar field problem, i.e. single unknown per node (DOF = 1)
					nodeValuePair.insert(std::pair<unsigned int, double>((dof * (dbNode - 1) + j), dirichletBValues[i]));
				}
			}
		}
	}

	// Impose values in final solution
	Eigen::MatrixXd result = ImposeDBC(solutionRed, nodeValuePair, numNodesDomain, dof);

	return result;
}


