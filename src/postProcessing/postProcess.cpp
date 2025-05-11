
#include "postProcess.hpp"



void GenerateVTK(std::string fileName, const std::vector<std::vector<double>> &NODE_COORD, const Element &meshElement, const Equation &equation, const Eigen::MatrixXd &Solution)
{
    // Writting VTK output file
    std::ofstream file_vtk;
    file_vtk.open(fileName);

    Eigen::IOFormat HeavyFmt(Eigen::FullPrecision);

    if (!file_vtk){
    	std::cerr<<"PostProcessing ERROR : Not able to write data to VTK File.\n Couldn't open file"<<std::endl;
    	exit (-403);
    }

    file_vtk << "# vtk DataFile Version 2.0" << "\n";
    //file_vtk << "output file at time "<< inputs.TotalTime <<std::endl;

    // Title of file
    file_vtk << "output test file" << "\n";

    // format for subsequent data
    file_vtk << "ASCII" << std::endl;

    // Dataset type
    file_vtk << "DATASET UNSTRUCTURED_GRID" << "\n";

    // Nodal coodinates
    file_vtk << "POINTS "<< meshElement.numNodes <<" double" << "\n";
    for (unsigned int i=0; i< meshElement.numNodes; i++){
	    file_vtk<<std::fixed << std::setprecision(8)<< NODE_COORD[i][0] << " " << NODE_COORD[i][1] << " " << NODE_COORD[i][2] << "\n";
    }

    int cellsType;

    // Linear triangle 2D element (3 node)
    if (meshElement.elemType == 2){
        file_vtk << "CELLS "<<meshElement.numElems<<" "<< meshElement.numElems*4 <<"\n";
        cellsType = 5;
    }
    // Linear quadratic 2D element (4 node)
    else if(meshElement.elemType == 3){
        file_vtk << "CELLS "<<meshElement.numElems<<" "<< meshElement.numElems*5 <<"\n";
        cellsType = 9;
        // cellsType = 8;
    }
    // Tetrahedron (4 node) 3D element
    else if (meshElement.elemType == 4){
        file_vtk << "CELLS "<<meshElement.numElems<<" "<< meshElement.numElems*5 <<"\n";
        cellsType = 10;
    }
    // Hexahedron 3D element (8 Nodes)
    else if (meshElement.elemType == 5){
        file_vtk << "CELLS "<<meshElement.numElems<<" "<< meshElement.numElems*9 <<"\n";
        cellsType = 12;
    }
    // Triangle second order 2D element (6 Nodes)
    else if (meshElement.elemType == 9){
        file_vtk << "CELLS "<<meshElement.numElems<<" "<< meshElement.numElems*7 <<"\n";
        cellsType = 22;
    }
    // Quadrangle second order 2D element (9 Nodes)
    else if (meshElement.elemType == 10){
        file_vtk << "CELLS "<<meshElement.numElems<<" "<< meshElement.numElems*10 <<"\n";
        cellsType = 23;
    }


    // Export all CELLS data [Element connectivity]
    for (unsigned int i = 0; i < meshElement.numElems; i++){
        file_vtk << meshElement.numNodesElem;
        for (int j = 0, nj = meshElement.ElemConnectivity[i].size(); j < nj; j++){
            file_vtk << " " << meshElement.ElemConnectivity[i][j] - 1;
        }
        file_vtk << "\n";
    }


    file_vtk << "CELL_TYPES "<<meshElement.numElems<<std::endl;
    file_vtk << Eigen::VectorXi::Ones(meshElement.numElems) * cellsType << "\n";


    file_vtk << "POINT_DATA "<< meshElement.numNodes<<"\n";

    // Thermal results
    if (equation.solverEq == 2)
        file_vtk << "SCALARS Temperature double 1"<<"\n";
    
    // Displacement results [2D]
    else if (equation.solverEq == 1 || equation.solverEq == 13){
        file_vtk << "SCALARS Displacement double 2" << "\n";
        file_vtk << "LOOKUP_TABLE default" << "\n";
        file_vtk << Solution.format(HeavyFmt) << "\n";

        file_vtk.close();
        return;
    }
    // Displacement results [3D]
    else if (equation.solverEq == 14){
        file_vtk << "SCALARS Displacement double 3" << "\n";
    }

    file_vtk << "LOOKUP_TABLE default" << "\n";
    file_vtk << Solution.format(HeavyFmt) <<std::endl;
    file_vtk.close();
}


/*
// Function to get Solution values at specific time intervals
void getTemporalResult(Eigen::MatrixXd Solution, double start, double end, readMesh mesh, std::string fileType, std::string unknown, SolverInp solverInp)
{

    // Convert all characters to upper to eliminate user input discrepancy
    for (int i = 0, n = fileType.length(); i<n; i++)
	{
		fileType[i] = toupper(fileType[i]);
	}

    // Character array to generate name of the file autonomously
    char outFileName[40];

    std::cout << "#### Writing solution to Results file ####" << std::endl;
    // std::cout << "Type of problem = " << inputs.problemType << std::endl;
    // std::cout << "Solution values for Mesh = " << mesh.fileName << std::endl;
    // std::cout << "Solution at each node is in ascending order of node numbers i.e. 1 to N.\n" << std::endl;

    // Columns counter for solution, i.e. get the column of solution for respective time step selected by user.
    int j = ceil((start - solverInp.StartTime)/solverInp.dt);

    // If need results for a single time step
    if (start == end)
    {
        if (fileType == "CSV"){
            snprintf(outFileName, sizeof(outFileName), "results/Result_%.3lf.csv", end);
            ExportCSV(Solution.col(j), outFileName);
        }
        else{
            snprintf(outFileName, sizeof(outFileName), "results/Result_%.3lf.txt", end);
            ExportTXT(Solution.col(j), outFileName, unknown);
        }
    }
    else
    {

        for (int i = 0, n = ceil((end-start)/solverInp.dt); i <= n; i++)
        {
            double time = start + solverInp.dt*i;

            if (fileType == "CSV"){
                snprintf(outFileName, sizeof(outFileName), "results/Result_%.3lf.csv", time);
                ExportCSV(Solution.col(j), outFileName);
            }
            else{
                snprintf(outFileName, sizeof(outFileName), "results/Result_%.3lf.txt", time);
                ExportTXT(Solution.col(j), outFileName, unknown);
            }
            j++;
        }
    }
}
*/


// Function to generate a CSV file of results
void ExportCSV(const Eigen::MatrixXd &Solution, std::string fileName)
{
    // File output stream object
    std::ofstream file;
    file.open(fileName);

    // Checking for file error, terminate if can't write out data to file
    if (!file){
        std::cerr<<"ERROR: Not able to write data to a CSV file."<<std::endl;
        std::cerr<<" Could not open file."<<std::endl;
        std::cerr<<" Filename: "<<fileName<<std::endl;
        std::cerr<<" Help:"<<std::endl;
        std::cerr<<"  Kindly try to verify that the filepath specified in the filename exists."<<std::endl;
        std::cerr<<"  If the filepath (i.e folder structure) specified in the filename above does not exist, then create it."<<std::endl;
        exit(1);
    }

    std::cout << "## Writing out data to a csv file = " << fileName <<std::endl;

    // const Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",","\n");
    const Eigen::IOFormat CSVFormat(10, Eigen::DontAlignCols, ",", "\n");
    file << Solution.format(CSVFormat);

    file.close();
}



// Write out solution to a .txt file
void ExportTXT(const Eigen::VectorXd &Solution,std::string outFileName, std::string name = "")
{
    std::ofstream file;
    file.open(outFileName);

    // Checking for file error, terminate if can't write out data to file
    if (!file){
        std::cerr<<"ERROR: Not able to write data to a TXT file."<<std::endl;
        std::cerr<<" Could not open file."<<std::endl;
        std::cerr<<" Filename: "<<outFileName<<std::endl;
        std::cerr<<" Help:"<<std::endl;
        std::cerr<<"  Kindly try to verify that the filepath specified in the filename exists."<<std::endl;
        std::cerr<<"  If the filepath (i.e folder structure) specified in the filename above does not exist, then create it."<<std::endl;
        exit(1);
    }

    std::cout << "## Writing out solution to a txt file = " << outFileName <<" ##\n";

    file << name << std::endl;
    file << Solution;

    file.close();
}






// Stress calculation post processing
Eigen::MatrixXd Stress_PostProcess(const Eigen::VectorXd &Displacements, std::vector<std::vector<unsigned int>> Domain_Connect, const std::vector<std::vector<double>> &NODE_COORD, const Eigen::MatrixXd &Emat, const readMesh &mesh)
{
    Eigen::MatrixXd Stress;

    ShapeFn2D BasisFn;
    BasisFn.getShapeFn(mesh.elemType);

    // Jacobian Matrix
    Eigen::Matrix2d Jacobian;

    // Element Connectivity vector
    std::vector<unsigned int> Te(Domain_Connect[0].size());

    for (int ie = 0, nie = Domain_Connect.size(); ie < nie; ie++)
    {
        Te = Domain_Connect[ie];

        // Element Node Coordinates
        Eigen::MatrixXd Element_NC(Te.size(), mesh.NODE_DIM);
        Element_NC.setZero();

        int ct = 0;
        for (unsigned int x : Te)
        {
            if (mesh.NODE_DIM == 3)
            {
                Element_NC.row(ct) << NODE_COORD[x-1][0], NODE_COORD[x-1][1], NODE_COORD[x-1][2];
            }
            else{
                Element_NC.row(ct) << NODE_COORD[x-1][0], NODE_COORD[x-1][1];
            }
            ct++;
        }

        for (int i = 0; i < BasisFn.NGP; i++)
        {
            Eigen::VectorXd Nx(BasisFn.NGP);
            Eigen::VectorXd Ny(BasisFn.NGP);

            // Calculate element Jacobian used to convert from natural (isoparametric) to global (Cartesian) coordinates.
            Jacobian = Jacobian2D(Nx, Ny, Element_NC, BasisFn, i);

            // Strain-Displacement matrix and Matrix of Shape function
            Eigen::MatrixXd B, Nmat;
        
            // 3 node, linear triangle element
            if (BasisFn.eType == 2)
            {
                B = Eigen::MatrixXd::Zero(3, 6);
                B.row(0) << Nx(0), 0, Nx(1), 0, Nx(2), 0;
                B.row(1) << 0, Ny(0), 0, Ny(1), 0, Ny(2);
                B.row(2) << Ny(0), Nx(0), Ny(1), Nx(1), Ny(2), Nx(2);

                Nmat = Eigen::MatrixXd::Zero(2, 6);
                Nmat.row(0) << BasisFn.N(i,0), 0, BasisFn.N(i,1), 0, BasisFn.N(i,2),0;
                Nmat.row(1) << 0, BasisFn.N(i,0), 0, BasisFn.N(i,1), 0, BasisFn.N(i,2);
            }
            // 4 node, linear quadrilateral element
            else if (BasisFn.eType == 3)
            {
                B = Eigen::MatrixXd::Zero(3, 8);
                B.row(0) << Nx(0), 0, Nx(1), 0, Nx(2), 0, Nx(3), 0;
                B.row(1) << 0, Ny(0), 0, Ny(1), 0, Ny(2), 0, Ny(3);
                B.row(2) << Ny(0), Nx(0), Ny(1), Nx(1), Ny(2), Nx(2), Ny(3), Nx(3);

                Nmat = Eigen::MatrixXd::Zero(2, 8);
                Nmat.row(0) << BasisFn.N(i,0), 0, BasisFn.N(i,1), 0, BasisFn.N(i,2), 0, BasisFn.N(i,3), 0;
                Nmat.row(1) << 0, BasisFn.N(i,0), 0, BasisFn.N(i,1), 0, BasisFn.N(i,2), 0, BasisFn.N(i,3);

            }
            // 6 node, second order triangle element
            else if (BasisFn.eType == 9)
            {
                B = Eigen::MatrixXd::Zero(3, 12);
                B.row(0) << Nx(0), 0, Nx(1), 0, Nx(2), 0, Nx(3), 0, Nx(4), 0, Nx(5), 0;
                B.row(1) << 0, Ny(0), 0, Ny(1), 0, Ny(2), 0, Ny(3), 0, Ny(4), 0, Ny(5);
                B.row(2) << Ny(0), Nx(0), Ny(1), Nx(1), Ny(2), Nx(2), Ny(3), Nx(3), Ny(4), Nx(4), Ny(5), Nx(5);

                Nmat = Eigen::MatrixXd::Zero(2, 12);
                Nmat.row(0) << BasisFn.N(i,0), 0, BasisFn.N(i,1), 0, BasisFn.N(i,2), 0, BasisFn.N(i,3), 0, BasisFn.N(i,4), 0, BasisFn.N(i,5), 0;
                Nmat.row(1) << 0, BasisFn.N(i,0), 0, BasisFn.N(i,1), 0, BasisFn.N(i,2), 0, BasisFn.N(i,3), 0, BasisFn.N(i,4), 0, BasisFn.N(i,5);
            }
            // Quadrilateral element, 2nd order {9 noded, quadratic interpolation}
            else if (BasisFn.eType == 10)
            {
                B = Eigen::MatrixXd::Zero(3, 18);
                B.row(0) << Nx(0), 0, Nx(1), 0, Nx(2), 0, Nx(3), 0, Nx(4), 0, Nx(5), 0, Nx(6), 0, Nx(7), 0, Nx(8), 0 ;
                B.row(1) << 0, Ny(0), 0, Ny(1), 0, Ny(2), 0, Ny(3), 0, Ny(4), 0, Ny(5), 0, Ny(6), 0, Ny(7), 0, Ny(8);
                B.row(2) << Ny(0), Nx(0), Ny(1), Nx(1), Ny(2), Nx(2), Ny(3), Nx(3), Ny(4), Nx(4), Ny(5), Nx(5), Ny(6), Nx(6), Ny(7), Nx(7), Ny(8), Nx(8);

                Nmat = Eigen::MatrixXd::Zero(2, 18);
                Nmat.row(0) << BasisFn.N(i,0), 0, BasisFn.N(i,1), 0, BasisFn.N(i,2), 0, BasisFn.N(i,3), 0, BasisFn.N(i,4), 0, BasisFn.N(i,5), 0, BasisFn.N(i,6), 0, BasisFn.N(i,7), 0, BasisFn.N(i,8), 0;
                Nmat.row(1) << 0, BasisFn.N(i,0), 0, BasisFn.N(i,1), 0, BasisFn.N(i,2), 0, BasisFn.N(i,3), 0, BasisFn.N(i,4), 0, BasisFn.N(i,5), 0, BasisFn.N(i,6), 0, BasisFn.N(i,7), 0, BasisFn.N(i,8);
            }
        }
    }

    return Stress;
}



