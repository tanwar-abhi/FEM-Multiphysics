
#include "Boundary_HT.hpp"
#include "structureBoundary.hpp"


// Heat load vector contribution from Heat source/sink boudary i.e. line heat source
// void LineHeatSource(Eigen::VectorXd &fGlobal, std::vector<std::vector<double>> NODE_COORD, readMesh mesh, ProbParameters_HT input, std::vector<std::vector<unsigned int>> T)
void LineHeatSource(Eigen::VectorXd &fGlobal, std::vector<std::vector<double>> NODE_COORD, readMesh mesh, std::vector<std::vector<unsigned int>> T, SolverInp solverInp, MaterialThermal materialDetails, NeumannBC boundary)
{

    // // Length of the Boundary/Surface Edge where heat generation source/sink exist
    // double LEdge;
    // // Calculate Edge length of the whole boundary
    // LEdge = sqrt(pow((NODE_COORD[T[T.size()-1][1]-1][0] - NODE_COORD[T[0][0]-1][0] ), 2)  +  pow((NODE_COORD[T[T.size()-1][1]-1][1] - NODE_COORD[T[0][0]-1][1] ), 2));


    // Number of element node
    int nen;

    if (mesh.BeType == 1){
        nen = 2;
    }
    else if (mesh.BeType == 8){
        nen = 3;
    }

    ShapeFn1D BasisFn1D;
    BasisFn1D.getShapeFn(mesh.BeType, nen);

    // Current element connectivity
    std::vector<unsigned int> Te(nen);

    // Vector Contribution due to Surface heat Conduction
    Eigen::VectorXd f_Q(nen);


    // Loop over each element of boundary
    for (int i = 0, n = T.size(); i < n; i++)
    {
        Te = T[i];
        f_Q.setZero();

        // Element node coordinates
        Eigen::MatrixXd Xe(Te.size(), mesh.NODE_DIM);
        Xe.setZero();

        int ct = 0;

        for (unsigned int x : Te)
        {
            if (mesh.NODE_DIM == 3)
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


        if (solverInp.coordinateSystem == "AXIS")
        {
            /*
            // Heat load vector contribution due to internal heat generation (RQ)
            if (value.domainSource && value.Q != 0){
                Eigen::VectorXd termF(nen);
                termF.row(0) << 2*Xe(0,0) + Xe(1,0) + Xe(2,0);
                termF.row(1) << Xe(0,0) + 2*Xe(1,0) + Xe(2,0);
                termF.row(2) << Xe(0,0) + Xe(1,0) + 2*Xe(2,0);

                flocal += value.Q * (area_elem/12.0) * termF;
            }
            */
        }
        else{
            // Loop over integration points of 1D element
            for (int j = 0; j < BasisFn1D.NGP; j++)
            {
                // Heat load vector contribution due to heat source/sink
                // In 1d element det(Jacobian) = le/2
                // f_Q += input.Q * BasisFn1D.N.row(j) * input.thickness * BasisFn1D.WGP[j] * le/2;
                f_Q += boundary.Q * BasisFn1D.N.row(j) * materialDetails.thickness * BasisFn1D.WGP[j] * le/2;
            }
        }

        for (int j = 0, kl = Te.size(); j < kl; j++)
        {
            // Adding element traction to global force vector
            fGlobal(Te[j]-1) += f_Q(j);
        }
    }
}




/*/ Heat load vector contribution due to radiant heat flow
void SurfaceRadiation(ShapeFn1D BasisFn1D, Eigen::VectorXd &fGlobal, std::vector<std::vector<double>> NODE_COORD, const float thickness, readMesh mesh, const ProblemParameters::HeatTransfer_inp input)
{
    // Boundary nodes connectivity of elements of surface where radiation heat flow rate is defined 
    // Connectivity of boundary element nodes
    std::vector<std::vector<unsigned int>> T;

    // Length of the Neuman Boundary Edge
    double LEdge;

    // Number of element node;{currently consider only linear elements}
    int nen = 2; 

    // Calculate Edge length of the boundary element
    LEdge = sqrt(pow((NODE_COORD[T[T.size()-1][1]-1][0] - NODE_COORD[T[0][0]-1][0] ), 2)  +  pow((NODE_COORD[T[T.size()-1][1]-1][1] - NODE_COORD[T[0][0]-1][1] ), 2));

    // if (BeType == 1){
    //     nen = 2;
    // }
    // else if (BeType == 8){
    //     nen = 3;
    // }

    // Current element connectivity
    std::vector<unsigned int> Te(nen);

    // Vector Contribution due to Surface convection
    Eigen::VectorXd R_r(nen); R_r.setZero();

    // Loop over each element of neuman boundary
    for (int i = 0, n = T.size(); i < n; i++)
    {
        Te = T[i];

        // Element node coordinates
        Eigen::MatrixXd Xe(Te.size(), mesh.NODE_DIM);
        Xe.setZero();
        int ct = 0;

        for (unsigned int x : Te)
        {
            if (mesh.NODE_DIM == 3)
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

        // Initializing traction force vector for element
        Eigen::VectorXd fe_Surf(nen);
        fe_Surf.setZero();


        if (input.problemType == "AXSYM")
        {
            // Integration for axisymmetric equations
        }
        else{
            // Loop over integration points of 1D element
            for (int j = 0; j < BasisFn1D.NGP; j++)
            {
                // Surface Radiation due to radiant heat flow over boundary
                R_r += input.ALPHA * input.qs * BasisFn1D.N.row(j) * le/2 * BasisFn1D.WGP[j];
            }
        }

        for (int j = 0, kl = Te.size(); j < kl; j++)
        {
            // Adding element traction to global force vector
            fGlobal(Te[j]-1) += R_r(j);
        }
    }
}



// Heat load vector contribution due to radiant heat flow
void SurfaceRadiation_HF(ShapeFn1D BasisFn1D, Eigen::VectorXd &fGlobal, std::vector<std::vector<double>> NODE_COORD, readMesh mesh, std::vector<std::vector<unsigned int>> T, SolverInp solverInp, NeumannBC boundary, MaterialThermal matDetails)
{
    // Boundary nodes connectivity of elements of surface where radiation heat flow rate is defined 
    // Connectivity of boundary element nodes
    // std::vector<std::vector<unsigned int>> T;

    // // Length of the Neuman Boundary Edge
    // double LEdge;
    // // Calculate Edge length of the boundary element
    // LEdge = sqrt(pow((NODE_COORD[T[T.size()-1][1]-1][0] - NODE_COORD[T[0][0]-1][0] ), 2)  +  pow((NODE_COORD[T[T.size()-1][1]-1][1] - NODE_COORD[T[0][0]-1][1] ), 2));

    // Number of element node;{currently consider only linear elements}
    int nen = 2; 


    if (mesh.BeType == 1){
        nen = 2;
    }
    else if (mesh.BeType == 8){
        nen = 3;
    }

    // Current element connectivity
    std::vector<unsigned int> Te(nen);

    // Vector Contribution due to Surface convection
    Eigen::VectorXd R_r(nen); 

    // Loop over each element of neuman boundary
    for (int i = 0, n = T.size(); i < n; i++)
    {
        Te = T[i];
        R_r.setZero();

        // Element node coordinates
        Eigen::MatrixXd Xe(Te.size(), mesh.NODE_DIM);
        Xe.setZero();
        int ct = 0;

        for (unsigned int x : Te)
        {
            if (mesh.NODE_DIM == 3)
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

        // Initializing traction force vector for element
        Eigen::VectorXd fe_Surf(nen);
        fe_Surf.setZero();


        // Outward normal vector
        std::vector<double> normal = OutwardNormal(Xe, le);

        if (input.problemType == "AXSYM")
        {
            Eigen::VectorXd term(Te.size());
            term.row(0) << 2*Xe(0,0) + Xe(1,0);
            term.row(1) << Xe(0,0) + 2*Xe(1,0);

            // Surface radiation contribution to heat load vector
            // R_r += (input.ALPHA * (input.qr[0]*normal[0] + input.qr[1]*normal[1]) * le/6.0) * term;
            R_r += (boundary.Radiation.ALPHA * (boundary.Radiation.qr[0]*normal[0] + boundary.Radiation.qr[1]*normal[1]) * le/6.0) * term;
        }
        else{
            // Loop over integration points of 1D element
            for (int j = 0; j < BasisFn1D.NGP; j++)
            {
                // Surface Radiation due to radiant heat flow over boundary
                R_r += boundary.Radiation.ALPHA * (boundary.Radiation.qr[0]*normal[0] + boundary.Radiation.qr[1]*normal[1])*le * BasisFn1D.N.row(j) * BasisFn1D.WGP[j] * le/2;
            }
        }

        for (int j = 0, kl = Te.size(); j < kl; j++)
        {
            // Adding element traction to global force vector
            fGlobal(Te[j]-1) += R_r(j);
        }
    }
}
*/



// Unit normal for each element
std::vector<double> outwardNormal(const Eigen::MatrixXd &Elem_Coord, double elemLength, int elementType)
{
    // Unit normal vector in for of [nx, ny] or [nx, ny, nz]
    std::vector<double> normalVector;

    // 2D domain, i.e. boundary has 1D line elements
    if (elementType == 1 || elementType == 8)
    {
        // double lengthB = sqrt( pow( (Elem_Coord(1,0)-Elem_Coord(0,0)), 2) );
        // double angle = asin(lengthB / elemLength);
        // normalVector.push_back(cos(angle));
        // normalVector.push_back(sin(angle));

        double length = sqrt(pow((Elem_Coord(1,0)-Elem_Coord(0,0)), 2) + pow((Elem_Coord(1,1)-Elem_Coord(0,1)), 2));
        normalVector.push_back((Elem_Coord(1,1)-Elem_Coord(0,1))/length);
        normalVector.push_back((Elem_Coord(0,0)-Elem_Coord(1,0))/length);
    }
    // 3D domain, with 2D triangle or quadrilateral elements at boundary
    else if (elementType == 2 || elementType == 9 || elementType == 3 || elementType == 10){
        Eigen::Vector3d v1 = Elem_Coord.row(1) - Elem_Coord.row(0); 
        Eigen::Vector3d v2 = Elem_Coord.row(2) - Elem_Coord.row(0);
        Eigen::Vector3d crossProduct = v1.cross(v2);
        crossProduct /= crossProduct.norm();
        for (int i = 0; i < 3; i++){
            normalVector.push_back(crossProduct(i));
        }
    }
    else{
        std::cerr<<"Error in OutwardNormal\n";
        std::cerr<<"Boundary element surface normal not defined \n";
        exit(-405);
    }

    return normalVector;
}






// ############### NEW improved implementations #############

// Heat load vector contribution from Surface conduction due to heat flux over surface area (S2)
void SurfaceConduction(Eigen::VectorXd &fGlobal, const std::vector<std::vector<double>> &NODE_COORD, const Element &meshElement, 
                        const SolverInp &solverInp, const NeumannBC &boundary, const MaterialThermal &material)
{

    // Number of element node;
    int nen = meshElement.numNodesElem; 

    // Current element connectivity
    std::vector<unsigned int> Te(nen);

    // Vector Contribution due to Surface heat Conduction
    Eigen::VectorXd Rq(nen);
    Eigen::VectorXd heatFlux(solverInp.dimension);
    // Traction/Pressure/Force components vector




    // Loop over each element of boundary
    for (int i = 0, n = meshElement.ElemConnectivity.size(); i < n; i++)
    {
        Te = meshElement.ElemConnectivity[i];
        Rq.setZero();

        // Element node coordinates
        Eigen::MatrixXd Xe(Te.size(), solverInp.dimension);
        Xe.setZero();

        int ct = 0;

        for (unsigned int x : Te)
        {
            if (solverInp.dimension == 3)
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

        // Unit Normal to suface
        std::vector<double> normal;
        normal = outwardNormal(Xe, le, meshElement.elemType);
        double heatFluxNormal = 0.0;
        if (boundary.boundaryType == "COMPONENTS"){
            for (int k = 0; k < solverInp.dimension; k++){
                    heatFluxNormal += boundary.values[k] * normal[k];
            }
        }
        else if (boundary.boundaryType == "NORMALTOBOUNDARY"){
            heatFluxNormal = boundary.value;
        }
        
        // Axisymetric problem type
        if (solverInp.coordinateSystem == "AXIS")
        {
            Eigen::VectorXd term(Te.size());
            term.row(0) << 2*Xe(0,0) + Xe(1,0);
            term.row(1) << Xe(0,0) + 2*Xe(1,0);

            // Integration for axisymmetric equations
            Rq += -((boundary.values[0]*normal[0] + boundary.values[1]*normal[1])*le * le/6.0) * term;
        }
        else{

            // Number of gauss points
            int ngp;
            if (solverInp.dimension == 2)
                ngp = meshElement.basisFn1D.NGP;
            else if (solverInp.dimension == 3){
                ngp = meshElement.basisFn2D.NGP;
            }
            // Loop over integration points of 1D element
            for (int j = 0; j < ngp; j++)
            {
                if (solverInp.dimension == 2){
                    // Surface Heating
                    // qs.n * Length of boundary * Shape function * Weights from Gauss Quadratures * det(Jacobian)
                    Rq += heatFluxNormal * meshElement.basisFn1D.N.row(j) * meshElement.basisFn1D.WGP[j] * material.thickness * le/2;
                }
                if (solverInp.dimension == 3){
                    // Surface Heating
                    // qs.n * Length of boundary * Shape function * Weights from Gauss Quadratures * det(Jacobian)
                    double detJ = temp_detJ(Xe);
                    Rq += heatFluxNormal * meshElement.basisFn2D.N.row(j) * meshElement.basisFn2D.WGP[j] * detJ;
                }

            }
        }

        for (int j = 0, kl = Te.size(); j < kl; j++)
        {
            // Adding element traction to global force vector
            fGlobal(Te[j]-1) += Rq(j);
        }
    }
}





// Heat load vector contribution from Surface convection over surface area (S3)
void SurfaceConvection(Eigen::VectorXd &fGlobal, const std::vector<std::vector<double>> &NODE_COORD, const Element &meshElement, 
                        const SolverInp &solverInp, const NeumannBC &boundary, const MaterialThermal &material)
{

    // Number of element node;
    int nen = meshElement.numNodesElem; 

    // Current element connectivity
    std::vector<unsigned int> Te(nen);

    // Vector Contribution due to Surface convection
    Eigen::VectorXd Rh(nen);


    // Loop over each element of convection boundary
    for (int i = 0, n = meshElement.ElemConnectivity.size(); i < n; i++)
    {
        // Connectivity of current element i.e. in loop
        Te = meshElement.ElemConnectivity[i];
        Rh.setZero();

        // Element node coordinates
        Eigen::MatrixXd Xe(Te.size(), solverInp.dimension);
        Xe.setZero();
        int ct = 0;

        for (unsigned int x : Te)
        {
            if (solverInp.dimension == 3)
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


        if (solverInp.coordinateSystem == "AXIS")
        {
            Eigen::VectorXd term(Te.size());
            term.row(0) << 2*Xe(0,0) + Xe(1,0);
            term.row(1) << Xe(0,0) + 2*Xe(1,0);

            // Surface convection contribution to heat load vector
            Rh += (boundary.H * boundary.ambientTemp * le/6.0) * term;
        }
        else{
            // Number of gauss points
            int ngp;
            if (solverInp.dimension == 2)
                ngp = meshElement.basisFn1D.NGP;
            else if (solverInp.dimension == 3){
                ngp = meshElement.basisFn2D.NGP;
            }
            // Loop over integration points of 1D element
            for (int j = 0; j < ngp; j++)
            {
                if (solverInp.dimension == 2){
                    // Surface Convection 
                    Rh += boundary.H * boundary.ambientTemp * meshElement.basisFn1D.N.row(j) 
                            * meshElement.basisFn1D.WGP[j] *  material.thickness * le/2;
                }
                if (solverInp.dimension == 3){
                    // Surface Heating
                    // qs.n * Length of boundary * Shape function * Weights from Gauss Quadratures * det(Jacobian)
                    double detJ = temp_detJ(Xe);
                    Rh += boundary.H * boundary.ambientTemp * meshElement.basisFn2D.N.row(j) * meshElement.basisFn2D.WGP[j] * detJ;
                }

            }
        }

        for (int j = 0, kl = Te.size(); j < kl; j++)
        {
            // Adding element traction to global force vector
            fGlobal(Te[j]-1) += Rh(j);
        }
    }
}





// Stiffness matrix contribution from convection boundary (S3)
void StiffnessConvection(std::vector<Eigen::Triplet<double>> &TrList, const std::vector<std::vector<double>> &NODE_COORD, 
                            const Element &meshElement, const SolverInp &solverInp, const NeumannBC &boundary, const MaterialThermal &material)
{

    // Number of element node;
    int nen = meshElement.numNodesElem;

    // Current element connectivity
    std::vector<unsigned int> Te(nen);

    // Contribution due to Surface convection
    Eigen::MatrixXd Kh(nen,nen);

    // Loop over each element of the convection boundary (S3)
    for (int i = 0, n = meshElement.ElemConnectivity.size(); i < n; i++)
    {
        // Geet current element connectivity of loop
        Te = meshElement.ElemConnectivity[i];
        Kh.setZero();

        // Element node coordinates
        Eigen::MatrixXd Xe(Te.size(), solverInp.dimension);
        Xe.setZero();

        int ct = 0;

        for (unsigned int x : Te)
        {
            if (solverInp.dimension == 3)
            {
                Xe.row(ct) << NODE_COORD[x-1][0], NODE_COORD[x-1][1], NODE_COORD[x-1][2];
            }
            else{
                Xe.row(ct) << NODE_COORD[x-1][0], NODE_COORD[x-1][1];
            }
            ct++;
        }

        // Length of element on the convection boundary
        double le;

        if (meshElement.elemType == 1){
            le = sqrt(pow((Xe(Xe.rows()-1, 0) - Xe(0,0)), 2) + pow((Xe(Xe.rows()-1, 1) - Xe(0,1)), 2));
        }
        else if (meshElement.elemType == 8)
        {
            le = sqrt(pow((Xe(1,0) - Xe(0,0)), 2) + pow((Xe(1,1) - Xe(0,1)), 2));
        }


        if (solverInp.coordinateSystem == "AXIS")
        {
            // Multiplication term for axisymmetric
            Eigen::MatrixXd term;

            if (meshElement.elemType == 1){
                term = Eigen::MatrixXd::Zero(2,2);
                term.row(0) << 3*Xe(0,0)+Xe(1,0) , Xe(0,0)+Xe(1,0);
                term.row(1) << Xe(0,0)+Xe(1,0) , Xe(0,0)+3*Xe(1,0);
            }
            else if (meshElement.elemType == 8){
                // term = Eigen::MatrixXd::Zero(3,3);
                std::cerr<<"\n # Higher order element not available in axisymmetric case # \n";
                exit(-405);
            }

            // Integration for axisymmetric equations
            Kh += boundary.H * (le/12.0) * term;
        }
        else{

            // Number of gauss points
            int ngp = 0;

            // 3D domain means 2D boundary
            if (solverInp.dimension == 3){
                ngp = meshElement.basisFn2D.NGP;
            }
            // 2D domain corresponds to 1D boundary
            else if (solverInp.dimension == 2){
                ngp = meshElement.basisFn1D.NGP;
            }

            // Loop over integration points of 1D element
            for (int j = 0; j < ngp; j++)
            {
                // Surface Convection
                if (solverInp.dimension == 2){
                    Kh += boundary.H * meshElement.basisFn1D.N.row(j).transpose() * meshElement.basisFn1D.N.row(j) * le / 2.0 
                            *  material.thickness * meshElement.basisFn1D.WGP[j];
                }
                else if (solverInp.dimension == 3){
                    double detJ = temp_detJ(Xe);
                    Kh += boundary.H * meshElement.basisFn2D.N.row(j).transpose() * 
                            meshElement.basisFn2D.N.row(j) * meshElement.basisFn2D.WGP[j] * detJ;
                }
            }
        }

        // std::cout<<"at neuman element ("<<i<<") stiffness convection =\n"<<Kh<<std::endl;

        // Add each convection boundary element contribution to global stiffness matrix
        for (int m = 0, lm = Kh.rows(); m < lm; m++){
            for (int n = 0, ln = Kh.cols(); n < ln; n++)
            {
                TrList.push_back(Tr(Te[m]-1, Te[n]-1, Kh(m,n)));
            }
        }
    }
}


