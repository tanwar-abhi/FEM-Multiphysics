
#include "structureBoundary.hpp"


/*/ Initialize vector of boundary force terms
void initializeBForce(Eigen::VectorXd &BForce, Parameters_LE inputs)
{
    if (inputs.bForce[0] == 0 && inputs.bForce[1] == 0)
    {
        BForce.setZero();
    }
    else{
        // BForce = Eigen::VectorXd::Ones(BForce.size());
        for (int i = 0, n = BForce.size()/2; i < n; i++)
        {
            if (i == 0){
                BForce(i) = inputs.bForce[i];
            }
            else{
                BForce(2*i) = inputs.bForce[1];
                BForce(2*i + 1) = inputs.bForce[0];
            }
        }
    }
}
*/



/* Imposing neuman boundary forces term in global force vector
void imposeTraction(ShapeFn1D shapeFunction1D, Eigen::VectorXd &fGlobal, std::vector<std::vector<double>> NODE_COORD,std::vector<std::vector<unsigned int>> T, readMesh mesh, SolverInp solverInpObj, BoundaryConditions boundaryTraction)
{

    // Length of the Neuman Boundary Edge (for calculating tractions)
    double LEdge;

    // Number of element node;
    int nen;

    // Traction force matrix
    // NOTE :: For only 2D is defined change accordingly for 3D
    Eigen::Vector2d trF;
    // trF(0) = tForce[0]; trF(1) = tForce[1];

    // Calculate Edge length of the neuman boundary
    LEdge = sqrt(pow((NODE_COORD[T[T.size()-1][1]-1][0] - NODE_COORD[T[0][0]-1][0] ), 2)  +  pow((NODE_COORD[T[T.size()-1][1]-1][1] - NODE_COORD[T[0][0]-1][1] ), 2));

    if (mesh.BeType == 1){
        nen = 2;
    }
    else if (mesh.BeType == 8){
        nen = 3;
    }

    // Current element connectivity
    std::vector<unsigned int> Te(nen);

    // Initializing traction force contribution to load for element
    Eigen::VectorXd fe_tr(nen*2);


    // Current Element node coordinates
    // Loop over each element of neuman boundary
    for (int i = 0, n = T.size(); i < n; i++)
    {
        Te = T[i];
        fe_tr.setZero();

        // Element node coordinates
        Eigen::MatrixXd Xe(Te.size(), solverInpObj.dimension);
        Xe.setZero();
        int ct = 0;

        for (unsigned int x : Te)
        {
            if (solverInpObj.dimension == 3)
            {
                Xe.row(ct) << NODE_COORD[x-1][0], NODE_COORD[x-1][1], NODE_COORD[x-1][2];
            }
            else{
                Xe.row(ct) << NODE_COORD[x-1][0], NODE_COORD[x-1][1];
            }
            ct++;
        }

        // Length of element
        double le;
        if (nen == 2){
            le = sqrt(pow((Xe(Xe.rows()-1, 0) - Xe(0,0)), 2) + pow((Xe(Xe.rows()-1, 1) - Xe(0,1)), 2));
        }
        else if (nen == 3)
        {
            le = sqrt(pow((Xe(1,0) - Xe(0,0)), 2) + pow((Xe(1,1) - Xe(0,1)), 2));
        }

        // Normal vector from surface
        std::vector<double> normal = OutwardNormal(Xe, le);

        // Traction force along y-axis thus angle of components would change
        if (input.tForce[0] == 0 && input.tForce[1] != 0)
        {
            trF(0) = input.tForce[0] * normal[1];
            trF(1) = input.tForce[1] * normal[0];
        }
        else{
            // Traction force components as per orientation of surface using normal vectors
            for (int ti = 0; ti < mesh.DIM; ti++)
            {
                trF(ti) = input.tForce[ti] * normal[ti];
            }
        }


        if (input.problemType == "AXSYM")
        {
            double a, b;
            a = (2*Xe(0,0) + Xe(Xe.rows()-1 ,0))/6;
            b = (Xe(0,0) + 2*Xe(Xe.rows()-1 ,0))/6;

            Eigen::VectorXd TrF_AxSym(nen*2);
            TrF_AxSym.row(0) << a*trF[0]; TrF_AxSym.row(1) << a*trF[1];
            TrF_AxSym.row(2) << b*trF[0]; TrF_AxSym.row(3) << b*trF[1];

            fe_tr += 2 * M_PI * le * TrF_AxSym;
        }
        else{
            // Loop over integration points of 1D element
            for (int j = 0; j < shapeFunction1D.NGP; j++)
            {
                Eigen::MatrixXd Nmat(2, nen*2);
                Nmat.setZero();
                if (nen == 2){
                    Nmat.row(0) << shapeFunction1D.N(j,0), 0, shapeFunction1D.N(j,1), 0;
                    Nmat.row(1) << 0, shapeFunction1D.N(j,0), 0, shapeFunction1D.N(j,1);
                }
                else if(nen == 3){
                    Nmat.row(0) << shapeFunction1D.N(j,0), 0, shapeFunction1D.N(j,1), 0, shapeFunction1D.N(j,2), 0;
                    Nmat.row(1) << 0, shapeFunction1D.N(j,0), 0, shapeFunction1D.N(j,1), 0, shapeFunction1D.N(j,2);
                }

                // Boundary element traction force calculation.
                fe_tr += input.thickness * (Nmat.transpose() * (trF/LEdge) * shapeFunction1D.WGP[j] * le/2);
                // fe_tr += input.thickness * (Nmat.transpose() * (trF/le) * shapeFunction1D.WGP[j] * le/2);
            }
        }

        // Local element numbers to global nodes numbers, as per 2dof
        std::vector<unsigned int> Te2dof(nen*2);
        ct = 0;
        for (int j = 0 , kl = nen; j < kl; j++)
        {            
            Te2dof[ct] = 2*(Te[j] - 1);
            Te2dof[ct+1] = 2*(Te[j] - 1) + 1;
            ct += 2;
        }

        for (int j = 0, kl = Te2dof.size(); j < kl; j++)
        {
            // Adding element traction to global force vector
            fGlobal(Te2dof[j]) += fe_tr(j);
        }   
    }
}
*/

//temporary function to calculate the detJ. Need to be implemented based on Jacobinan based method
double temp_detJ(const Eigen::MatrixXd &Xe)
{
	Eigen::Vector3d AB = Xe.row(1) - Xe.row(0);
	Eigen::Vector3d AC = Xe.row(2) - Xe.row(0);
	Eigen::Vector3d AD = Xe.row(3) - Xe.row(0);
	double surfaceArea = 0.5* (AC.cross(AB).norm() + AC.cross(AD).norm());
    return surfaceArea/4;
}

//Function to find surface area
double getSurfaceArea(const Element &meshElement, const readMesh &mesh, const SolverInp &solverInpObj, 
                        const LinearElasticMaterial &material)
{
    double surfaceArea = 0.0;
    // Loop over each element of traction boundary
    for (int i = 0, n = meshElement.numElems; i < n; i++)
    {
        // Current element connectivity (within elemental loop)
        std::vector<unsigned int> Te = meshElement.ElemConnectivity[i];

        // Element node coordinates
        Eigen::MatrixXd Xe(Te.size(), solverInpObj.dimension);
        Xe.setZero();
        int ct = 0;
        for (unsigned int x : Te)
        {
            if (solverInpObj.dimension == 3){
                Xe.row(ct) << mesh.Node_Coord[x-1][0], mesh.Node_Coord[x-1][1], mesh.Node_Coord[x-1][2];
            }
            else{
                Xe.row(ct) << mesh.Node_Coord[x-1][0], mesh.Node_Coord[x-1][1];
            }
            ct++;
        }
        if (solverInpObj.dimension == 2){
            double LEdge;
            LEdge = sqrt(pow((Xe(1,0) - Xe(0,0)), 2)  +  pow((Xe(1,1) - Xe(0,1)), 2));

            surfaceArea += LEdge*material.thickness;
        }
        else if (solverInpObj.dimension == 3){

		    Eigen::Vector3d AB = Xe.row(1) - Xe.row(0);
		    Eigen::Vector3d AC = Xe.row(2) - Xe.row(0);
		    Eigen::Vector3d AD = Xe.row(3) - Xe.row(0);
		    surfaceArea += 0.5* (AC.cross(AB).norm() + AC.cross(AD).norm());
        }
    }

    return surfaceArea;
}



// Imposing neuman boundary forces term in global force vector
void imposeTraction(Eigen::VectorXd &fGlobal, const Element &meshElement, const readMesh &mesh, const SolverInp &solverInpObj, const NeumannBC &boundary, const LinearElasticMaterial &material, int eqn)
{
    double trFNormal = 0.0;
    Eigen::VectorXd trF(solverInpObj.dimension);

    // Traction/Pressure/Force components vector
    if (boundary.variable == "TRACTION" || boundary.variable == "PRESSURE" || boundary.variable == "FORCE") 
    {
        if (boundary.boundaryType == "COMPONENTS"){
            for (int i = 0; i < solverInpObj.dimension; i++){
                trF(i) = boundary.values[i];
            }
        }
        else if (boundary.boundaryType == "NORMALTOBOUNDARY"){
            trFNormal = boundary.value;
        }
    }

    double surfaceArea = getSurfaceArea(meshElement, mesh, solverInpObj, material);
    // Initializing traction force contribution to load for element
    Eigen::VectorXd fe_tr(meshElement.numNodesElem * solverInpObj.equations[eqn].DOF);


    // Loop over each element of traction boundary
    for (int i = 0, n = meshElement.numElems; i < n; i++)
    {
        // Current element connectivity (within elemental loop)
        std::vector<unsigned int> Te = meshElement.ElemConnectivity[i];

        // Reset element traction vector
        fe_tr.setZero();

        // Element node coordinates
        Eigen::MatrixXd Xe(Te.size(), solverInpObj.dimension);
        Xe.setZero();

        int ct = 0;
        for (unsigned int x : Te)
        {
            if (solverInpObj.dimension == 3){
                Xe.row(ct) << mesh.Node_Coord[x-1][0], mesh.Node_Coord[x-1][1], mesh.Node_Coord[x-1][2];
            }
            else{
                Xe.row(ct) << mesh.Node_Coord[x-1][0], mesh.Node_Coord[x-1][1];
            }
            ct++;
        }

        // Length of element
        double le = 0.0;
        // Linear element 2 node (2D)
        if (meshElement.elemType == 1){
            le = sqrt(pow((Xe(Xe.rows()-1, 0) - Xe(0,0)), 2) + pow((Xe(Xe.rows()-1, 1) - Xe(0,1)), 2));
        }
        // 3-node second order line element (2D)
        else if (meshElement.elemType == 8){
            le = sqrt(pow((Xe(1,0) - Xe(0,0)), 2) + pow((Xe(1,1) - Xe(0,1)), 2));
        }

        if (boundary.boundaryType == "NORMALTOBOUNDARY"){        
            // Normal vector from surface
            std::vector<double> normal = outwardNormal(Xe, le, meshElement.elemType);
            for (int ti = 0; ti < solverInpObj.dimension; ti++){
                // trF(ti) = boundary.neumann[indexTraction].values[ti] * normal[ti];
                trF(ti) =  normal[ti]*trFNormal;
            } 
        }

        if (solverInpObj.coordinateSystem == "AXSYM")
        {
            double a, b;
            a = (2*Xe(0,0) + Xe(Xe.rows()-1 ,0))/6;
            b = (Xe(0,0) + 2*Xe(Xe.rows()-1 ,0))/6;

            Eigen::VectorXd TrF_AxSym(meshElement.numNodesElem*2);
            TrF_AxSym.row(0) << a*trF[0]; TrF_AxSym.row(1) << a*trF[1];
            TrF_AxSym.row(2) << b*trF[0]; TrF_AxSym.row(3) << b*trF[1];

            fe_tr += 2 * M_PI * le * TrF_AxSym;
        }
        else{

            // Number of gauss points
            int ngp;
            if (solverInpObj.dimension == 2)
                ngp = meshElement.shapeFunction1D.NGP;
            else if (solverInpObj.dimension == 3){
                ngp = meshElement.shapeFunction2D.NGP;
            }

            Eigen::MatrixXd Jacobian(2,2); Jacobian.setZero();

            // Eigen::MatrixXd Nmat;
            // if (solverInpObj.dimension == 2)
            //     Nmat = Eigen::MatrixXd::Zero(2, meshElement.numNodesElem * 2);
            // else if (solverInpObj.dimension == 3){
            //     Nmat = Eigen::MatrixXd::Zero(3, meshElement.numNodesElem * solverInpObj.equations[eqn].DOF);
            // }

            // Loop over integration points of boundary element
            for (int j = 0; j < ngp; j++)
            {
                // Eigen::MatrixXd Nmat(2, meshElement.numNodesElem * 2);
                Eigen::MatrixXd Nmat(solverInpObj.dimension, meshElement.numNodesElem * solverInpObj.equations[eqn].DOF);


                // 2Node line element (1D) on boundary
                if (meshElement.elemType == 1)
                {
                    Nmat.row(0) << meshElement.shapeFunction1D.N(j,0), 0, meshElement.shapeFunction1D.N(j,1), 0;
                    Nmat.row(1) << 0, meshElement.shapeFunction1D.N(j,0), 0, meshElement.shapeFunction1D.N(j,1);
                }
                // 3-node second order line (1D) on boundary
                else if (meshElement.elemType == 8 )
                {
                    Nmat.row(0) << meshElement.shapeFunction1D.N(j,0), 0, meshElement.shapeFunction1D.N(j,1), 0, meshElement.shapeFunction1D.N(j,2), 0;
                    Nmat.row(1) << 0, meshElement.shapeFunction1D.N(j,0), 0, meshElement.shapeFunction1D.N(j,1), 0, meshElement.shapeFunction1D.N(j,2);
                }
                // 3-node triangle (2D) on boundary
                else if (meshElement.elemType == 2)
                {
                    // Nmat.row(0) << meshElement.shapeFunction2D.N(j,0), 0, meshElement.shapeFunction2D.N(j,1), 0, meshElement.shapeFunction2D.N(j,2), 0;
                    // Nmat.row(1) << 0, meshElement.shapeFunction2D.N(j,0), 0, meshElement.shapeFunction2D.N(j,1), 0, meshElement.shapeFunction2D.N(j,2);

                    Nmat.row(0) << meshElement.shapeFunction2D.N(j,0), 0, 0, meshElement.shapeFunction2D.N(j,1), 0, 0, meshElement.shapeFunction2D.N(j,2), 0, 0;
                    Nmat.row(1) << 0, meshElement.shapeFunction2D.N(j,0), 0, 0, meshElement.shapeFunction2D.N(j,1), 0, 0, meshElement.shapeFunction2D.N(j,2), 0;
                    Nmat.row(2) << 0, 0, meshElement.shapeFunction2D.N(j,0), 0, 0, meshElement.shapeFunction2D.N(j,1), 0, 0, meshElement.shapeFunction2D.N(j,2);
                }
                // 4-node rectangle (2D) on boundary
                else if (meshElement.elemType == 3){
                    // Nmat.row(0) << meshElement.shapeFunction2D.N(j,0), 0, meshElement.shapeFunction2D.N(j,1), 0, meshElement.shapeFunction2D.N(j,2), 0, meshElement.shapeFunction2D.N(j,3), 0;
                    // Nmat.row(1) << 0, meshElement.shapeFunction2D.N(j,0), 0, meshElement.shapeFunction2D.N(j,1), 0, meshElement.shapeFunction2D.N(j,2), 0, meshElement.shapeFunction2D.N(j,3);

                    Nmat.row(0) << meshElement.shapeFunction2D.N(j,0), 0, 0, meshElement.shapeFunction2D.N(j,1), 0, 0, meshElement.shapeFunction2D.N(j,2), 0, 0, meshElement.shapeFunction2D.N(j,3), 0, 0;
                    Nmat.row(1) << 0, meshElement.shapeFunction2D.N(j,0), 0, 0, meshElement.shapeFunction2D.N(j,1), 0, 0, meshElement.shapeFunction2D.N(j,2), 0, 0, meshElement.shapeFunction2D.N(j,3), 0;
                    Nmat.row(2) << 0, 0, meshElement.shapeFunction2D.N(j,0), 0, 0, meshElement.shapeFunction2D.N(j,1), 0, 0, meshElement.shapeFunction2D.N(j,2), 0, 0, meshElement.shapeFunction2D.N(j,3);
                }
                // Boundary element traction force calculation.
                if (solverInpObj.dimension == 2){
                    if (boundary.variable == "TRACTION" || boundary.variable == "PRESSURE")
                        fe_tr += material.thickness * (Nmat.transpose() * trF * meshElement.shapeFunction1D.WGP[j] * le / 2);
                    else if (boundary.variable == "FORCE")
                        fe_tr += material.thickness * (Nmat.transpose() * (trF / surfaceArea) * meshElement.shapeFunction1D.WGP[j] * le / 2);
                }
                else if (solverInpObj.dimension == 3){
                    // Derivative of shape functions (N) w.r.t global axis i.e.. X,Y,Z
                    Eigen::VectorXd Nx(ngp);
                    Eigen::VectorXd Ny(ngp);

                    //Jacobian = Jacobian2D2(Nx, Ny, Xe, meshElement.shapeFunction2D, j);

                    //double detJ = Jacobian.determinant();
                    //detJ = abs(detJ);
                    // Temporary method for determinant. Need to verify for imperfect quadrilateral 
                    double detJ = temp_detJ(Xe);
                    // fe_tr += (Nmat.transpose() * (trF/le) * shapeFunction1D.WGP[j] * le/2);
                    if (boundary.variable == "TRACTION" || boundary.variable == "PRESSURE")
                        fe_tr += (Nmat.transpose() * trF * meshElement.shapeFunction2D.WGP[j] * detJ);
                    else if (boundary.variable == "FORCE"){
                        fe_tr += (Nmat.transpose() * (trF / surfaceArea) * meshElement.shapeFunction2D.WGP[j] * detJ);
                    }
                }
            }
        }

        // Local element node numbers to global nodes numbers, as per DOF for problem
        std::vector<unsigned int> Te2dof(meshElement.numNodesElem * solverInpObj.equations[eqn].DOF);
        ct = 0;
        for (int j = 0, kl = meshElement.numNodesElem; j < kl; j++){
            for(int k = 0; k < solverInpObj.equations[eqn].DOF; k++){
                Te2dof[ct] = solverInpObj.equations[eqn].DOF * (Te[j] - 1) + k;
                ct++;
            }
        }

        // Adding element traction contribution to global force vector
        for (int j = 0, kl = Te2dof.size(); j < kl; j++){
            fGlobal(Te2dof[j]) += fe_tr(j);
        }
    }
}


