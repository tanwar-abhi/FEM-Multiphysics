
#include "Time_Discretization.hpp"
// #include "../Libraries/unsupported/Eigen/IterativeSolvers"
#include <Eigen/IterativeLinearSolvers>
// #include <cmath>

// Initial boundary conditions used for the first time step (i.e. t=0)
void initialConditions(Eigen::VectorXd &Temperature, std::vector<unsigned int> T0Nodes, double T0)
{

    // Imposing initial condtions of the nodes
    for (auto it = T0Nodes.begin(); it != T0Nodes.end(); it++)
    {
        Temperature.row(*it-1) << T0;
    }
}


/*
// Time discretizaiton scheme (Theta method aka Generalized Trapezoidal Rule)
// void TimeDiscretization_ThetaHT(Eigen::MatrixXd &Temp_sol, SpMatDouble C, SpMatDouble K, Eigen::VectorXd fglobal,std::vector<std::vector<double>> NODE_COORD, ProbParameters_HT inputs, std::vector<std::vector<unsigned int>> T0nodes, readMesh mesh, std::set<unsigned int> DB_Nodes)
void TimeDiscretization_ThetaHT(Eigen::MatrixXd &Temp_sol, SpMatDouble C, SpMatDouble K, Eigen::VectorXd fglobal, std::vector<std::vector<unsigned int>> EConnectivity, std::vector<std::vector<double>> NODE_COORD, std::vector<std::vector<unsigned int>> T0nodes, readMesh mesh, std::set<unsigned int> DB_Nodes, SolverInp solverInp, NeumannBC boundary, MaterialThermal material)
{

    // User defined dirichlet boundary conditions
    DirichletBC dirichlet;
    dirichlet.readBoundaryInputs();

    // user defined initial boundary conditions
    InitialBC intialBoundary;
    intialBoundary.readBoundaryInputs();


    // Create 1D shape function element to used by the boundary elements
    ShapeFn1D BasisFn1D;
    BasisFn1D.getShapeFn(mesh.BeType);

    // Vector of solution i.e. temperature at time step (n) {previous time step}
    Eigen::VectorXd temp_n(fglobal.size());

    // Temperature at time step (n+1) {current time step}
    Eigen::VectorXd temp_n1(temp_n.size());

    // This is only for linear problems; {LHS term}
    // for non linear equations element K,C,f must be evaluated at each time step
    SpMatDouble K_bar;
    K_bar = (1.0/solverInp.dt) * C + solverInp.THETA * K;

    std::vector<Eigen::VectorXd> moveToRHS;
    if (mesh.DBtag != -1 ){
        // Reduce stiffness Matrix and force vector to impose dirichlet boundary condition
        moveToRHS = EliminateLHSMatrix(K_bar, DB_Nodes, 1);

        // ####### Check this step??? is it only should be done once ##
        for (int i = 0, n = moveToRHS.size(); i < n; i++)
        {
            fglobal -= moveToRHS[i] * dirichlet.NodeValue;
        }
    }

    // Total iteration counter for each step
    unsigned int iteration;

    // Solve linear system of equation using ConjugateGradient solver
    Eigen::ConjugateGradient<SpMatDouble> solver;

    // Heat load vector (RHS term) at previous time step (n)
    Eigen::VectorXd F_n = fglobal;

    // Heat load vector (RHS term) at current time step (n+1)
    Eigen::VectorXd F_n1 = fglobal;

    // Running temporal loop for total no. of time steps
    for (int i = 0, n = ceil(solverInp.TotalTime/solverInp.dt) ; i < n; i++)
    {
        // First time step, thus apply initial condition (i.e. t = 0)
        if (i == 0)
        {
            temp_n.setZero(); temp_n1.setZero();
            initialConditions(temp_n, T0nodes, intialBoundary.IntialVale);
            temp_n = ImposeDBC(temp_n, DB_Nodes, dirichlet.NodeValue, mesh.NNodes);
            Temp_sol.col(i) << temp_n;
            // continue;
        }
        else{
            temp_n1.setZero();
            temp_n = Temp_sol.col(i-1);

            // Setting F_n = fglobal at previous timestep (n)
            // F_n = F_n1;
            // F_n1 = fglobal;
        }

        if (!boundary.HeatSource_Sink.domainSource && boundary.HeatSource_Sink.Q != 0)
        {
            std::vector<std::vector<unsigned int>> heatSourceNodes;
            heatSourceNodes = ReadBoundaryNodes(mesh.partsInp.part[0].meshFileName, "source");

            // LineHeatSource(fglobal, NODE_COORD, mesh, inputs, heatSourceNodes);
            LineHeatSource(fglobal, NODE_COORD, mesh, heatSourceNodes, solverInp, boundary, material);
        }

        // If Surface Heat flux defined (at conduction boundary) then calculate contribution to heat vector
        if (boundary.Conduction.qs[0] != 0 || boundary.Conduction.qs[1] != 0)
        {
            // Read nodes on the boundary where surface heating is defined
            std::vector<std::vector<unsigned int>> conductionBNodes;
            conductionBNodes = ReadBoundaryNodes(mesh.partsInp.part[0].meshFileName, "cond");

            // Heat load vector contribution from Surface heating over surface area (S2)
            // SurfaceConduction(BasisFn1D, F_n1, NODE_COORD, mesh, inputs, conductionBNodes);
            SurfaceConduction(BasisFn1D, F_n1, NODE_COORD, mesh, conductionBNodes, solverInp, boundary);
        }

        // Surface convection term contribution
        if (boundary.Convection.H != 0)
        {
            // Read nodes of elements on the boundary where convection is defnied
            std::vector<std::vector<unsigned int>> convectionBNodes;
            convectionBNodes = ReadBoundaryNodes(mesh.partsInp.part[0].meshFileName, "conv");

            // Also update the boundary name at line 504 at Matrix_Assemble
            // Heat load vector contribution from convection boundary over surface (S3)
            // SurfaceConvection(BasisFn1D, F_n1, NODE_COORD, mesh, inputs, convectionBNodes);
            SurfaceConvection(BasisFn1D, F_n1, NODE_COORD, mesh, convectionBNodes, solverInp, boundary);
        }


        // Surface radiation contribution due to radiant heat flow
        // if (inputs.qr[0] != 0 || inputs.qr[1] != 0)
        // {
        //     // Connectivity of Nodes on the radiation boundary where radiant heat flow acts
        //     std::vector<std::vector<unsigned int>> radiationBNodes;
        //     radiationBNodes = ReadBoundaryNodes(mesh.fileName, "PhysicalName of the radiation bounadry as per user");

        //     // Heat vector contribution due to radiation heat flow
        //     SurfaceRadiation_HF(BasisFn1D, F_n1, NODE_COORD, mesh, inputs, radiationBNodes);
        // }


        // Equivalent stiffness matrix and load vector from the recurance realtion
        // SpMatDouble K_bar(K.rows(), K.cols());
        Eigen::VectorXd f_bar(K.rows());

        iteration = 0;

        // intial error to start the while loop
        double error = 1;

        std::cout<<"#### Time step = " << i <<"  ####"<<std::endl;

        while (error > solverInp.eps){

            if (iteration != 0) {
                temp_n = temp_n1;
            }

            // Heat Load vector term {total RHS term}
            // Since we don't have any non linear term at the boundary thus heat load vector is same for current and previous time step {i.e. fn+1 = fn}
            // Re-check the term for non-linear problems, i.e. last 2 terms will not be fglobal but would be f_n and f_n1 i.e. heat load at current timestep and previous time step
            // f_bar = ((1.0/inputs.dt)*C - (1-inputs.THETA)*K)*temp_n + (1-inputs.THETA)*F_n + inputs.THETA*F_n1;
            f_bar = ((1.0/solverInp.dt)*C - (1-solverInp.THETA)*K)*temp_n + (1-solverInp.THETA)*F_n + solverInp.THETA*F_n1;


            if (mesh.DBtag != -1 ){
                for (int i = 0, n = moveToRHS.size(); i < n; i++)
                {
                    f_bar -= moveToRHS[i]* dirichlet.NodeValue;
                }
                // Reduce stiffness Matrix and force vector to impose dirichlet boundary condition
                // ReduceMatrix_1DOF(K_bar, f_bar, DB_Nodes);
                EliminateRHS(f_bar, DB_Nodes, 1);
            }

            // Solve linear system of equation using ConjugateGradient solver
            solver.compute(K_bar);

            temp_n1 = solver.solve(f_bar);

            // Impose i.e. replace the values at nodes of dirichlet boundary 
            if (mesh.DBtag != -1){
                temp_n1 = ImposeDBC(temp_n1, DB_Nodes, dirichlet.NodeValue, mesh.NNodes);
            }

            // Error vector for calculating l2 norm
            Eigen::VectorXd er = temp_n1 - temp_n;
            error = er.squaredNorm();

            if (iteration > 10000){
                std::cerr << "Convergance not reached within 10000 iterations at time step(i) = "<<i<<std::endl;
                std::cerr << "\n i.e. at Total time = "<< solverInp.StartTime+(solverInp.dt*i)<<"\nCheck the Time Discretization function"<<std::endl;
                exit(10);
            }
            iteration++;
            std::cout<<"Iteration = "<<iteration<<std::endl;

        }

        // if (mesh.DBtag != -1){
        //     temp_n1 = ImposeDBC(temp_n1, DB_Nodes, inputs.DBV, mesh.NNodes);
        // }

        // Solution converged then save values for current time step in temperature solution matrix
        Temp_sol.col(i) << temp_n1;

        // ###################### CHECK ####################################
        char outFileName[40];
	char vtkOutFileName[40];
        //snprintf(outFileName, sizeof(outFileName), "test_%.3lf.csv", solverInp.StartTime+(solverInp.dt*i));
	//snprintf(vtkOutFileName, sizeof(vtkOutFileName), "test_%.3lf.vtk", solverInp.StartTime+(solverInp.dt*i));
        //ExportCSV(Temp_sol.col(i),outFileName, "Temperature");
	//GenerateVTK(vtkOutFileName, NODE_COORD, EConnectivity, mesh, Temp_sol.col(i));

        // std::cout << "Solution converged fat time step = " << inputs.StartTime + (inputs.dt * i);
    }

    std::cout<<"Solution converged, results obtained for all time steps.\n";
}
*/





/*
    if (solverInpObj.transient)
        // Global conductance {mass} matrix
        SpMatDouble C(mesh.NNodes, mesh.NNodes);

        // Add a for loop based on the equation class obj
        std::cout << "Started Transient Simulations "<< solverInpObj.nEquations << std::endl;
        for (int eqn=0; eqn < solverInpObj.nEquations; eqn++)
        {
            std::cout << "Equation type is " << solverInpObj.equations[eqn].solverEq << std::endl;
            if (solverInpObj.equations[eqn].solverEq==2)
            {
                // std::cout << " Started Solving Heat Transfer " << std::endl;

                // Assemble the global stiffness matrix and Force vector for Linear Elasticity problem
                int elemTagId = solverInpObj.equations[eqn].elemTagId;

                assembleMatrixHT(K, C, f, mesh.Node_Coord, meshElements[elemTagId].ElemConnectivity, meshElements[elemTagId].basisFn2D, meshElements, elemTagId, solverInpObj, materialThermal, boundary);
            }
	    }
        // std::cout<<"F vector = \n"<<f<<std::endl;

        // Total iteration counter for each step
        unsigned int iteration;

        // Selecting ConjugateGradient solver
        Eigen::ConjugateGradient<SpMatDouble> solver;

        // Heat load vector (RHS term) at previous time step (n)
        Eigen::VectorXd F_n = f;

        // Heat load vector (RHS term) at current time step (n+1)
        Eigen::VectorXd F_n1 = f;

        // Vector of unknowns at previous time step (n) {temporary variable; thermal->temperature, solid->displacements}
        Eigen::VectorXd temp_n(f.size());

        // Vector of unknowns at current time step (n+1) {temporary variable; thermal->temperature, solid->displacements}
        Eigen::VectorXd temp_n1(temp_n.size());

        // This is only for linear problems; {LHS term}
        // for non linear equations element K,C,f must be evaluated at each time step
        SpMatDouble K_bar;
        K_bar = (1.0 / solverInpObj.dt) * C + solverInpObj.THETA * K;

        // Running temporal loop for total no. of time steps
        for (int i = 0, n = ceil(solverInpObj.TotalTime/solverInpObj.dt) ; i < n; i++)
        {
            // First time step, thus apply initial condition (i.e. t = 0)
            if (i == 0)
            {
                temp_n.setZero(); temp_n1.setZero();
                for (int index=0; index<boundary.nIBC; index++){
                    if (boundary.initial[index].type == 1){
                        int elemTagId = boundary.initial[index].elemTagId;
                        initialConditions(temp_n, meshElements[elemTagId].Nodes, boundary.initial[index].value);
                    }
                }
                // continue;
                F_n.setZero();
            }
            else{
                temp_n1.setZero();

                // Change this to get from the temporary variable and avoid saving whole time solution
                // Setting solution at previous time step for next iteration
                temp_n = Solution.col(i-1);

                // Setting F_n = fglobal at previous timestep (n)
                F_n = F_n1;
            }


            for (int index = 0; index < boundary.nNBC; index++)
            {

                // If Surface Heat flux defined (at conduction boundary) then calculate contribution to heat vector
                if (boundary.neuman[index].variable == "HEATFLUX")
                {
                    if (boundary.neuman[index].values[0] != 0 || boundary.neuman[index].values[1] != 0)
                    {
                        // Read nodes on the boundary where surface heating is defined
                        int elementTag = boundary.neuman[index].elemTagId;

                        // Heat load vector contribution from Surface heating over surface area (S2)
                        SurfaceConduction(meshElements[elementTag].basisFn1D, F_n1, mesh.Node_Coord, meshElements[elementTag], meshElements[elementTag].ElemConnectivity, solverInpObj, boundary.neuman[index]);
                    }
                }

                else if (boundary.neuman[index].variable == "CONVECTIVEHEATTRANSFER")
                {
                    // Surface convection term contribution
                    if (boundary.neuman[index].H != 0)
                    {
                        // Read nodes of elements on the boundary where convection is defnied
                        int elementTag = boundary.neuman[index].elemTagId;

                        // Heat load vector contribution from convection boundary over surface (S3)
                        SurfaceConvection(meshElements[elementTag].basisFn1D, F_n1, mesh.Node_Coord, meshElements[elementTag], meshElements[elementTag].ElemConnectivity, solverInpObj, boundary.neuman[index]);
                    }
                }

                // else if (boundary.neuman[index].domainSource == false && boundary.neuman[index].variable == "HEATSOURCESINK")
                // {
                //     if(boundary.neuman[index].Q != 0)
                //     {
                //         std::vector<std::vector<unsigned int>> heatSourceNodes = meshElements[boundary.neuman[index].elemTagId].ElemConnectivity;
                //         // std::vector<std::vector<double>> NODE_COORD = meshElements[boundary.neuman[index].elemTagId].Nodes;
                //         // heatSourceNodes = ReadBoundaryNodes(mesh.fileName, "source");

                //         LineHeatSource(f, mesh.Node_Coord, mesh, heatSourceNodes, solverInpObj, material_Thermal[0], boundary.neuman[index]);
                //     }
                // }

            }


            // Equivalent stiffness matrix and load vector from the recurance realtion
            // SpMatDouble K_bar(K.rows(), K.cols());
            Eigen::VectorXd f_bar(K_bar.rows()); f_bar.setZero();

            iteration = 0;

            // intial error to start the while loop
            double error = 3.14152;

            std::cout << "#### Time step = " << i << "  ####" << std::endl;

            std::cout<<"F_n = \n"<<F_n<<std::endl;
            std::cout<<"F_n1 = \n"<<F_n1<<std::endl;
            std::cout<<"temp_n = \n"<<temp_n<<std::endl;
            std::cout<<"temp_n1 = \n"<<temp_n1<<std::endl;

            while (error > solverInpObj.eps)
            {

                if (iteration != 0)
                {
                    temp_n = temp_n1;
                }

                // Heat Load vector term {total RHS term}

                // Pardha -- This step will be very challenging tomorrow when we want to make our code in to sparse.
                f_bar = ((1.0 / solverInpObj.dt) * C - (1 - solverInpObj.THETA) * K) * temp_n + (1 - solverInpObj.THETA) * F_n + solverInpObj.THETA * F_n1;

                std::cout<<"f_bar {time = "<<i<<",iteration "<<iteration<<"} = \n"<<f_bar<<std::endl;

                // Adding contribution to force vector
                if (iteration == 0)
                {
                    f_bar += f;
                }

                // F_n1 = f_bar;

                // ########################### change as per new logics #################.
                if (boundary.nDBC != 0)
                {
                    /*
                    // Vector containing element nodes of all the elements of the dirichlet boundary
                    std::vector<std::vector<unsigned int>> allDBnodes;

                    // Vector contaiing all values of the subsequent boundary nodes
                    std::vector<double> allDBvalues;

                    for (int index=0; index<boundary.nDBC; i++)
                    {
                        int elementTadId = boundary.dirichlet[index].elemTagId;
                        
                        allDBnodes.push_back(meshElements[elementTadId].Nodes);
                        allDBvalues.push_back(boundary.dirichlet[index].value);
                    }
                    EliminateMultiDBC(allDBnodes, allDBvalues, K_bar, f_bar, 1);

                }


                // Solve linear system of equation using ConjugateGradient solver
                solver.compute(K_bar);

                temp_n1 = solver.solve(f_bar);


                // Impose i.e. replace the values of nodes within solution for dirichlet boundary
                if (boundary.nDBC != 0)
                {
                    for (int i=0; i<boundary.nDBC; i++)
                    {
                        temp_n1 = ImposeDBC(temp_n1, meshElements[boundary.dirichlet[i].elemTagId].Nodes, boundary.dirichlet[i].value, mesh.NNodes);
                    }
                }

                // Error vector for calculating l2 norm
                std::cout<<"temp_n = \n"<<temp_n<<std::endl;
                std::cout<<"temp_n1 = \n"<<temp_n1<<std::endl;

                Eigen::VectorXd er = temp_n1 - temp_n;
                error = er.squaredNorm();

                if (iteration > 10000)
                {
                    std::cerr << "Convergance not reached within 10000 iterations at time step(i) = " << i << std::endl;
                    std::cerr << "\n i.e. at Total time = " << solverInpObj.StartTime + (solverInpObj.dt * i) << "\nCheck the Time Discretization function" << std::endl;
                    exit(10);
                }
                iteration++;
                std::cout << "Iteration = " << iteration << std::endl;
            }

            Solution.col(i) << temp_n1;

            // ###################### CHECK ####################################
            // char outFileName[40];
            // char vtkOutFileName[40];
            // snprintf(outFileName, sizeof(outFileName), "test_%.3lf.csv", solverInpObj.StartTime + (solverInpObj.dt * i));
            // snprintf(vtkOutFileName, sizeof(vtkOutFileName), "test_%.3lf.vtk", solverInpObj.StartTime + (solverInpObj.dt * i));
            // ExportCSV(Solution.col(i), outFileName, "Temperature");
            // GenerateVTK(vtkOutFileName, mesh.Node_Coord, meshElements[solverInpObj.equations[0].elemTagId].ElemConnectivity, mesh, Solution.col(i));

            // std::cout << "Solution converged fat time step = " << inputs.StartTime + (inputs.dt * i);
        }

        // Starting and ending point of the time step to get temporal data from Transient results
        double posStart, posEnd;
        std::cout << "Enter inital time to get temporal solution = ";
        std::cin >> posStart;

        std::cout << "Enter final time to get temporal solution = ";
        std::cin >> posEnd;

        getTemporalResult(Solution, posStart, posEnd, mesh, "csv", "Temperature", solverInpObj);

*/
