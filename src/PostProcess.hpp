// Post Process functions defined here to call upon after solver i.e. post calculation of results

#pragma once
#include <fstream>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include "Matrix_Assemble.hpp"


// Function to generate the vtk files that would be read by ParaView for post processing.
void GenerateVTK(std::string fileName, const std::vector<std::vector<double>> &NODE_COORD, const Element &meshElement, const Equation &equation, const Eigen::MatrixXd &Solution);


/* Function to obtain specific results from whole temporal solution
Solution = Matrix containing nodal values at each time step
start = Starting position of the time step from where temporal soultion is desired
end = Ending position of the time step i.e. upto which time the solutioin is desired
meshObject = object of the readMesh class containing mesh data
userInputs = onject or problem Parameters containing all constants and user inputs for the problem
fileeType = Type of file to store the solution, currently supports only CSV and TXT files
unknown = Name of the unknown (i.e. Displacement, Temperature, pressure etc.) to be added to solution file first row

void getTemporalResult(Eigen::MatrixXd Solution, double start, double end, readMesh meshObject, std::string fileType, std::string unknown, SolverInp solverInp);
*/

// Function to generate a CSV file of any variable
// This function can also be used in transient results by calling this function for each time step
void ExportCSV(const Eigen::MatrixXd &variable, std::string fileName);



/* Function to write out any variable to a .txt file
Solution = Vector containing nodal solution
fileName = Name of the results file with extension in string format i.e. fileNAME.txt
name = The name of the unknown i.e. name of data in csv (eg Displacement, Temperature, pressure etc.) in String format
*/
void ExportTXT(const Eigen::VectorXd &variable,std::string fileName, std::string name);



// Overload function to export coordinates and connectivity data from code to csv file
template<typename T>
void ExportCSV(std::vector<std::vector<T>> variable, std::string fileName)
{
    std::ofstream file;
    file.open(fileName);

    if (!file)
    {
        std::cerr<<"ExportCSV error :: ";
        std::cerr<<"Unable to open file to write csv data to file\n";
        exit(-405);
    }

    for (int i = 0, ni = variable.size(); i < ni; i++)
    {
        for (int j = 0, nj = variable[i].size(); j < nj; j++)
        {
            if (j == variable[i].size()-1)
            {
                file << variable[i][j] << std::endl;
            }
            else{
                file << variable[i][j] << ",";
            }
        }
    }

    file.close();
}

