
#include "Matrix_Assemble.hpp"

#include <Eigen/QR>
#include <stdexcept>


// To calculate centroid Point for Axisymmetric elements
std::vector<std::vector<double>> CentroidPoint(double &r_Bar, double &z_Bar, const Eigen::MatrixXd &Xe, int nen, int elemType );



// Calculates the Jacobian matrix and derivative of shape Function W.R.T x and y at current Gauss Point
Eigen::Matrix2d Jacobian2D(Eigen::VectorXd &Nx, Eigen::VectorXd &Ny, const Eigen::MatrixXd &Xe, const ShapeFunction2D &shapeFunction, const int itr)
{
    // Jacobian Matrix
    Eigen::Matrix2d J;
    J.setZero();

    J.row(0) << shapeFunction.Nxi.row(itr) * Xe.col(0) , shapeFunction.Nxi.row(itr) * Xe.col(1);
    J.row(1) << shapeFunction.Neta.row(itr) * Xe.col(0) , shapeFunction.Neta.row(itr) * Xe.col(1);

    // Matrix containing derivate of N w.r.t Natural coordinates (i.e. xi and eta)
    Eigen::MatrixXd NNatural(2, shapeFunction.Nxi.cols());
    NNatural.row(0) << shapeFunction.Nxi.row(itr);
    NNatural.row(1) << shapeFunction.Neta.row(itr);

    // Calculate derivative of Shape function w.r.t global coodinates (i.e. x, y)
    Eigen::MatrixXd result = J.completeOrthogonalDecomposition().solve(NNatural);

    Nx = result.row(0);
    Ny =  result.row(1);

    return J;
}




// Calculates the Jacobian matrix and derivative of shape Function of 3D element W.R.T x,y and z at current Gauss Point
Eigen::Matrix3d Jacobian3D(Eigen::VectorXd &Nx, Eigen::VectorXd &Ny, Eigen::VectorXd &Nz, const Eigen::MatrixXd &Xe, const ShapeFunction3D &shapeFunction, const int itr)
{
    // Jacobian Matrix
    Eigen::Matrix3d J;
    J.setZero();

    J.row(0) << shapeFunction.Nxi.row(itr) * Xe.col(0), shapeFunction.Nxi.row(itr) * Xe.col(1), shapeFunction.Nxi.row(itr) * Xe.col(2);
    J.row(1) << shapeFunction.Neta.row(itr) * Xe.col(0), shapeFunction.Neta.row(itr) * Xe.col(1), shapeFunction.Neta.row(itr) * Xe.col(2);
    J.row(2) << shapeFunction.Nzeta.row(itr) * Xe.col(0), shapeFunction.Nzeta.row(itr) * Xe.col(1), shapeFunction.Nzeta.row(itr) * Xe.col(2);

    // Matrix containing derivate of N w.r.t Natural coordinates (i.e. Xi, Eta, Zeta)
    Eigen::MatrixXd NNatural(3, shapeFunction.Nxi.cols());
    NNatural.row(0) << shapeFunction.Nxi.row(itr);
    NNatural.row(1) << shapeFunction.Neta.row(itr);
    NNatural.row(2) << shapeFunction.Nzeta.row(itr);

    // Calculate derivative of Shape function w.r.t global coodinates (i.e. x, y)
    Eigen::MatrixXd result = J.completeOrthogonalDecomposition().solve(NNatural);

    Nx = result.row(0);
    Ny = result.row(1);
    Nz = result.row(2);

    return J;
}



/* Function to calculate the Centroidal point of the element to be used by axysymmetric element
r_Bar = Mean radial distance from axis of symmetry
z_Bar = Mean vertical distance from axis of symmetry
Xe = Nodal Coodinates of the element
nen = No. of element nodes
elemType = Element type as per gmsh
*/
std::vector<std::vector<double>> CentroidPoint(double &r_Bar, double &z_Bar, const Eigen::MatrixXd &Xe, int nen, int elemType )
{

    std::vector<double> alpha, beta, gama;

    for (int i = 0; i<nen; i++)
    {
        r_Bar += Xe(i,0);
        z_Bar += Xe(i,1);

        if (elemType == 2){

            if (i+2 < nen){
                alpha.push_back( Xe(i+1,0)*Xe(i+2,1) - Xe(i+1,1)*Xe(i+2,0) );
                beta.push_back( Xe(i+1,1) - Xe(i+2,1) );
                gama.push_back( Xe(i+2,0) - Xe(i+1,0) );
            }
            else if (i+1 > nen-1)
            {
                alpha.push_back( Xe(i-2,0)*Xe(i-1,1) - Xe(i-2,1)*Xe(i-1,0) );
                beta.push_back( Xe(i-2,1) - Xe(i-1,1) );
                gama.push_back( Xe(i-1,0) - Xe(i-2,0) );
            }
            else{
                alpha.push_back( Xe(i+1,0)*Xe(0,1) - Xe(i+1,1)*Xe(0,0) );
                beta.push_back( Xe(i+1,1) - Xe(0,1) );
                gama.push_back( Xe(0,0) - Xe(i+1,0) );
            }
        }
        
        // Quadrilateral element
        // else if (elemType == 3)
        // {
            // area_elem = sqrt(pow(Xe(1,0)-Xe(0,0),2)+pow(Xe(1,1)-Xe(0,1),2)) * sqrt(pow(Xe(2,0)-Xe(1,0),2)+pow(Xe(2,1)-Xe(1,1),2));
        // }
    }

    // Centroidal point of element (r_Bar, z_Bar)
    r_Bar /= nen;
    z_Bar /= nen;

    std::vector<std::vector<double>> results;
    results.push_back(alpha);
    results.push_back(beta);
    results.push_back(gama);

    return results;
}





// Element Matrix for Heat transfer problem (Axisymetric and Cartesian)
void ElementMatrixHT(Eigen::MatrixXd &Klocal, Eigen::MatrixXd &Clocal, Eigen::VectorXd &flocal, const Eigen::MatrixXd &Xe, const Element &meshElement, const SolverInp &solverInp, const MaterialThermal &material, const BoundaryConditions &boundary)
{

    // Initialize the local stiffness matrix and source (RHS) vector to zero
    Klocal.setZero();
    Clocal.setZero();
    flocal.setZero();

    // The index at which heat source boundary details are in boundary array
    int indexHS = -404;
    for (int i = 0; i < boundary.nNBC; i++)
    {
        // Heat source term is similar to body force thus required by stiffness matrix function
        if (boundary.neumann[i].variable == "HEATSOURCESINK")
        {
            if (boundary.neumann[i].domainSource){
                indexHS = i;
                break;
            }
        }
    }

    // Jacobian Matrix
    Eigen::MatrixXd Jacobian;

    // Area of element
    double area_elem;

    // Centroidal points (used for axisymetrix element)
    double r_Bar, z_Bar;

    // variable for axisymmetric element 
    std::vector<double> alpha, beta, gama; 

    // Calculate element area and centroid (centroid only for axisymetric element)
    // Linear 3 nodes Triangle element
    if (meshElement.elemType == 2)
    {
        area_elem = (1.0/2.0) * (Xe(0,0)*(Xe(1,1)-Xe(2,1)) + Xe(1,0)*(Xe(2,1)-Xe(0,1)) + Xe(2,0)*(Xe(0,1)-Xe(1,1)));
        area_elem = sqrt(pow(area_elem,2));

        if (solverInp.coordinateSystem == "AXIS"){
            std::vector<std::vector<double>> ABG;

            // Calculate centroidal point for each element, alpha, beta and gama {ABG} variables for shape function
            ABG = CentroidPoint(r_Bar, z_Bar, Xe, meshElement.numNodesElem, meshElement.elemType);
            alpha = ABG[0]; beta = ABG[1]; gama = ABG[2];
        }
    }
    // Quadrilateral (Linear 4 node element)
    else if (meshElement.elemType == 3){
        area_elem = sqrt(pow(Xe(1,0)-Xe(0,0),2)+pow(Xe(1,1)-Xe(0,1),2)) * sqrt(pow(Xe(2,0)-Xe(1,0),2)+pow(Xe(2,1)-Xe(1,1),2));

        if (solverInp.coordinateSystem == "AXIS"){
            CentroidPoint(r_Bar, z_Bar, Xe, meshElement.numNodesElem, meshElement.elemType);
        }
    }

    if (solverInp.coordinateSystem == "AXIS"){
        // Triangle element (linear, 3 nodes)
        if (meshElement.elemType == 2)
        {
            // Derivative of shape functions W.R.T r and z
            std::vector<double> Nr, Nz;

            for (int i = 0; i < meshElement.numNodesElem; i++)
            {
                // N[i] = (1.0/(2*area_elem))*(alpha[i] + beta[i]*r_Bar + gama[i]*z_Bar);
                // N[i] = (1.0/(2*area_elem))*(alpha[i] + beta[i]*Xe(i,0) + gama[i]*Xe(i,1));
                Nr.push_back(beta[i] / (2 * area_elem));
                Nz.push_back(gama[i] / (2 * area_elem));
            }

            // Capacitance (Mass) Matrix for Transient problems
            if (solverInp.isTransient)
            {
                Eigen::MatrixXd term(meshElement.numNodesElem, meshElement.numNodesElem);
                term.row(0) << 6 * Xe(0, 0) + 2 * Xe(1, 0) + 2 * Xe(2, 0), 2 * Xe(0, 0) + 2 * Xe(1, 0) + Xe(2, 0), 2 * Xe(0, 0) + Xe(1, 0) + 2 * Xe(2, 0);
                term.row(1) << 2 * Xe(0, 0) + 2 * Xe(1, 0) + Xe(2, 0), 2 * Xe(0, 0) + 6 * Xe(1, 0) + 2 * Xe(2, 0), Xe(0, 0) + 2 * Xe(1, 0) + 2 * Xe(2, 0);
                term.row(2) << 2 * Xe(0, 0) + Xe(1, 0) + 2 * Xe(2, 0), Xe(0, 0) + 2 * Xe(1, 0) + 2 * Xe(2, 0), 2 * Xe(0, 0) + 2 * Xe(1, 0) + 6 * Xe(2, 0);

                Clocal += (material.RHO * material.spHeat * area_elem / 60.0) * term;
            }

            // Temperature gradient interpolation matrix
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2, meshElement.numNodesElem);
            B.row(0) << Nr[0], Nr[1], Nr[2];
            B.row(1) << Nz[0], Nz[1], Nz[2];

            // Thermal conductivity tensor for isotropic/anisotropic material
            Eigen::MatrixXd K_Tensor = Material_HT(B.rows(), material.k);

            // Conduciton matrix contribution to local stiffness matrix
            Klocal += r_Bar * area_elem * B.transpose() * K_Tensor * B;

            // If problem contains heat source/sink term only then add contributions to thermal force
            if (indexHS != -404)
            {
                // Heat load vector contribution due to internal heat generation (RQ)
                if (boundary.neumann[indexHS].domainSource && boundary.neumann[indexHS].Q != 0)
                {
                    Eigen::VectorXd termF(meshElement.numNodesElem);
                    termF.row(0) << 2 * Xe(0, 0) + Xe(1, 0) + Xe(2, 0);
                    termF.row(1) << Xe(0, 0) + 2 * Xe(1, 0) + Xe(2, 0);
                    termF.row(2) << Xe(0, 0) + Xe(1, 0) + 2 * Xe(2, 0);

                    flocal += boundary.neumann[indexHS].Q * (area_elem / 12.0) * termF;
                }
            }
        }
        // Higer element
        // else if (meshElement.elemType == 3){

        // }
    }
    // For cartesian coordinate system
    else
    {
        // Number of gauss points
        int ngp;
        if (solverInp.dimension == 2){
            ngp = meshElement.shapeFunction2D.NGP;
        }
        else if (solverInp.dimension == 3){
            ngp = meshElement.shapeFunction3D.NGP;
        }

        // Loop over gauss points (integration points)
        for (int i = 0; i < ngp; i++)
        {
            // Derivative of shape function N w.r.t "x", "y" and "z"
            Eigen::VectorXd Nx(ngp);
            Eigen::VectorXd Ny(ngp);
            Eigen::VectorXd Nz(ngp);

            // Calculate element Jacobian to convert from natural (isoparametric) to global (Cartesian) coordinates.
            if (solverInp.dimension == 2){
                Jacobian = Jacobian2D(Nx, Ny, Xe, meshElement.shapeFunction2D, i);
            }
            else if (solverInp.dimension == 3){
                Jacobian = Jacobian3D(Nx, Ny, Nz, Xe, meshElement.shapeFunction3D, i);
            }


            // Temperature gradient interpolation matrix
            Eigen::MatrixXd B(solverInp.dimension, meshElement.numNodesElem);
            B.setZero();


            // Linear Triangle element (3 nodes)
            if (meshElement.elemType == 2)
            {
                // meshElement.numNodesElem = 3;
                // B = Eigen::MatrixXd::Zero(2, 3);
                B.row(0) << Nx(0), Nx(1), Nx(2);
                B.row(1) << Ny(0), Ny(1), Ny(2);
            }
            // Linear Quadrilateral element (4 nodes)
            else if (meshElement.elemType == 3)
            {
                // meshElement.numNodesElem = 4;
                // B = Eigen::MatrixXd::Zero(2, 4);
                B.row(0) << Nx(0), Nx(1), Nx(2), Nx(3);
                B.row(1) << Ny(0), Ny(1), Ny(2), Ny(3);
            }
            // Triangle element (6 nodes) Quadratic interpolation
            else if (meshElement.elemType == 9)
            {
                // meshElement.numNodesElem = 6;
                // B = Eigen::MatrixXd::Zero(2, 6);
                B.row(0) << Nx(0), Nx(1), Nx(2), Nx(3), Nx(4), Nx(5);
                B.row(1) << Ny(0), Ny(1), Ny(2), Ny(3), Ny(4), Ny(5);
            }
            // Quadrilateral elemnt (9 nodes) Quadratic interpolation
            else if (meshElement.elemType == 10)
            {
                // meshElement.numNodesElem = 9;
                // B = Eigen::MatrixXd::Zero(2, 9);
                B.row(0) << Nx(0), Nx(1), Nx(2), Nx(3), Nx(4), Nx(5), Nx(6), Nx(7), Nx(8);
                B.row(1) << Ny(0), Ny(1), Ny(2), Ny(3), Ny(4), Ny(5), Ny(6), Ny(7), Ny(8);
            }
            // 8 Noded Hexahedral Linear element {3D}
            else if (meshElement.elemType == 5)
            {
                B.row(0) << Nx(0), Nx(1), Nx(2), Nx(3), Nx(4), Nx(5), Nx(6), Nx(7);
                B.row(1) << Ny(0), Ny(1), Ny(2), Ny(3), Ny(4), Ny(5), Ny(6), Ny(7);
                B.row(2) << Nz(0), Nz(1), Nz(2), Nz(3), Nz(4), Nz(5), Nz(6), Nz(7);
            }

            // Determinant of Jacobian
            double detJ = Jacobian.determinant();

            // Absolute of determinant (To eliminate negative determinant due to orientation of elements)
            // ################ To be changed later ##################
            detJ = abs(detJ);

            if (solverInp.isTransient)
            {
                // Consistent Mass (Capacitance) Matrix
                if (solverInp.massMatrixType)
                {
                    if (solverInp.dimension == 2)
                        Clocal += material.RHO * material.spHeat * meshElement.shapeFunction2D.N.row(i).transpose() * meshElement.shapeFunction2D.N.row(i) * material.thickness * meshElement.shapeFunction2D.WGP[i] * detJ;
                    else if (solverInp.dimension == 3){
                        Clocal += material.RHO * material.spHeat * meshElement.shapeFunction3D.N.row(i).transpose() * meshElement.shapeFunction3D.N.row(i) * material.thickness * meshElement.shapeFunction3D.WGP[i] * detJ;
                    }
                }
                // Lumped Capacitance Matrix Selected
                else
                {
                    Clocal += ((material.RHO * material.spHeat * material.thickness * area_elem) / 3.0) * Eigen::MatrixXd::Identity(meshElement.numNodesElem, meshElement.numNodesElem);
                }
            }

            // Thermal conductivity tensor for isotropic/anisotropic material
            Eigen::MatrixXd K_Tensor = Material_HT(B.rows(), material.k);
            if (solverInp.dimension == 2)
                Klocal += B.transpose() * K_Tensor * B * material.thickness * meshElement.shapeFunction2D.WGP[i] * detJ;
            else if (solverInp.dimension == 3){
                Klocal += B.transpose() * K_Tensor * B * meshElement.shapeFunction3D.WGP[i] * detJ;
            }


            // If problem contains heat source/sink term only then add contributions to thermal force vector
            if (indexHS != -404)
            {
                if (boundary.neumann[indexHS].domainSource && boundary.neumann[indexHS].Q != 0)
                {
                    // Internal heat generation contribution. (Rq)
                    if (solverInp.dimension == 2)
                        flocal += boundary.neumann[indexHS].Q * meshElement.shapeFunction2D.N.row(i) * material.thickness * meshElement.shapeFunction2D.WGP[i] * detJ;
                    else if (solverInp.dimension == 3){
                        flocal += boundary.neumann[indexHS].Q * meshElement.shapeFunction3D.N.row(i) * meshElement.shapeFunction3D.WGP[i] * detJ;
                    }
                }
            }
        }
    }
}





// ############################ Generalized expression for Linear elastic equation (2D + 3D)
// Using DOF variable defined in the equations class


/* Function to calculate element stiffness matrix and Force vector for Linear Elasticity
Input Arguments ::
    Klocal = Element stiffness Matrix
    flocal = Element force vector
    Xe = Element Node Coordintes, 
    Dmat = Material matrix based on hooke's law,
    meshElement = object of element containing domain elemenent connectivities and shape function
    solverInp = Object of SolverInp class containing sovler specific information
    material = object containing material properties
    boundary = Object containing all boundary conditions defined
*/
void ElementMatrixLE(Eigen::MatrixXd &K_local, Eigen::VectorXd &f_local, const Eigen::MatrixXd &Xe, const Eigen::MatrixXd &Dmat, const Element &meshElement, const SolverInp &solverInp, const LinearElasticMaterial &material, const BoundaryConditions &boundary)
{
    // Initialize all elements of element matrix/vector to zero
    K_local.setZero();
    f_local.setZero();

    // Index at which body force is defined
    int bodyForceIndex = -404;
    for (int i = 0; i < boundary.nDBC; i++){
        if (boundary.dirichlet[i].variable == "BODYFORCE")
        {
            bodyForceIndex = i;
            break;
        }
    }


    // Area of element
    double area_elem;

    // Axisymmetric Problem
    if (solverInp.coordinateSystem == "AXIS")
    {
        // Centroid point of the element (mean radial and axial distance of element)
        double r_Bar = 0.0;
        double z_Bar = 0.0;

        std::vector<std::vector<double>> ABG;

        // Calculate centroidal point for each element, alpha, beta and gama {ABG} variables for shape function
        ABG = CentroidPoint(r_Bar, z_Bar, Xe, meshElement.numNodesElem, meshElement.elemType);
        std::vector<double> alpha, beta, gama;
        alpha = ABG[0]; beta = ABG[1]; gama = ABG[2];

        // Shape function for Axisymmetric element
        std::vector<double> N(meshElement.numNodesElem);

        // Strain Displacement Matrix and Shape Function matrix
        Eigen::MatrixXd B, Nmat;

        // Linear Triangular axisymmetric element (2D)
        if (meshElement.elemType == 2){
            for (int i = 0; i < meshElement.numNodesElem; i++)
            {
                N[i] = (1.0/(2*area_elem))*(alpha[i] + beta[i]*r_Bar + gama[i]*z_Bar);
            }

            B = Eigen::MatrixXd::Zero(4, 6);
            B.row(0) << beta[0] , 0 , beta[1] , 0 , beta[2] , 0;
            B.row(1) << ((alpha[0]/r_Bar)+beta[0]+(gama[0]*z_Bar/r_Bar)) , 0 , ((alpha[1]/r_Bar)+beta[1]+(gama[1]*z_Bar/r_Bar)) , 0 , ((alpha[2]/r_Bar)+beta[2]+(gama[2]*z_Bar/r_Bar)) , 0;
            B.row(2) << 0 , gama[0] , 0 , gama[1] , 0 , gama[2];
            B.row(3) << gama[0], beta[0] , gama[1] , beta[1] , gama[2] , beta[2];
            B *= 1.0/(2*area_elem);

            Nmat = Eigen::MatrixXd::Zero(2,6);
            Nmat.row(0) << N[0] , 0 , N[1] , 0 , N[2] , 0;
            Nmat.row(1) << 0 , N[0] , 0 , N[1] , 0 , N[2];
        }

        // Linear Quadrilateral element (4 nodes) {2D}
        else if (meshElement.elemType == 3){
            B = Eigen::MatrixXd::Zero(4, 8);

            Nmat = Eigen::Matrix2d::Zero(2, 8);
            Nmat.row(0) << N[0], 0, N[1], 0, N[2], 0, N[3], 0;
            Nmat.row(1) << 0, N[0], 0, N[1], 0, N[2], 0, N[3];
        }


        // Jacobian
        Eigen::MatrixXd Jacob(2,2);
        Jacob.row(0) << Xe(0,0)-Xe(2,0) , Xe(0,1)-Xe(2,1);
        Jacob.row(1) << Xe(1,0)-Xe(2,0) , Xe(1,1)-Xe(2,1); 

        double detJ = (Xe(1,1)-Xe(2,1))*(Xe(0,0)-Xe(2,0)) - (Xe(0,1)-Xe(2,1))*(Xe(1,0)-Xe(2,0));
        detJ = abs(detJ);

        // Calculate local element stiffness matrix
        K_local += 2*M_PI * r_Bar * area_elem * B.transpose() * Dmat * B;


        // ############## Elemeent body forces to be moved or not to be moved ##############
        /*
        // If  Distributed body forces are defined(i.e not zero), only then calculate contribution to force vector
        if (!ElementBF.isZero())
        {
            // For a moving part Rb_bar would be non zero and Zb is body force per unit volume due to gravity
            double Rb_bar, Zb;
            if (ElementBF(1) != 0){
                Zb = -inputs.RHO*9.8*area_elem;
            }
            else{
                Zb = ElementBF(1);
            }

            Rb_bar = pow(inputs.omega,2)*inputs.RHO*r_Bar;
            
            f_local.row(0) << Rb_bar; f_local.row(1) << Zb;
            f_local.row(2) << Rb_bar; f_local.row(3) << Zb;
            f_local.row(4) << Rb_bar; f_local.row(5) << Zb;
            f_local *= (2*M_PI*r_Bar*area_elem/3.0);
        }
        */
    }

    // For problems in Cartesian coordinate
    else{

        // Jacobian Matrix
        Eigen::MatrixXd Jacobian;
        int numGaussPts;

        // Initialize variables as per 3D problem
        if (solverInp.dimension == 3)
        {
            numGaussPts = meshElement.shapeFunction3D.NGP;
            Jacobian = Eigen::MatrixXd::Zero(3,3);
        }
        // Initialize variables as per 2D problem
        else if (solverInp.dimension == 2)
        {
            numGaussPts = meshElement.shapeFunction2D.NGP;
            Jacobian = Eigen::MatrixXd::Zero(2,2);
        }

        // Loop over Gauss points (computation of integrals on the current element)
        for (int i = 0; i < numGaussPts; i++)
        {
            // Derivative of shape functions (N) w.r.t global axis i.e.. X,Y,Z
            // Nx = Derivative w.r.t X; Ny = Derivtive w.r.t Y; Nz = Derivative w.r.t Z;
            Eigen::VectorXd Nx(numGaussPts);
            Eigen::VectorXd Ny(numGaussPts);
            Eigen::VectorXd Nz(numGaussPts);

            // Calculate element Jacobian used to convert from natural (isoparametric) to global (Cartesian) coordinates.
            if (solverInp.dimension == 2)
                Jacobian = Jacobian2D(Nx, Ny, Xe, meshElement.shapeFunction2D, i);
            else if (solverInp.dimension == 3){
                Jacobian = Jacobian3D(Nx, Ny, Nz, Xe, meshElement.shapeFunction3D, i);
            }


            // B -> Strain-Displacement matrix ; Nmat -> Matrix of Shape function
            Eigen::MatrixXd B, Nmat;


            // 3 node, linear triangle element {2D}
            if (meshElement.elemType == 2)
            {
                B = Eigen::MatrixXd::Zero(3, 6);
                B.row(0) << Nx(0), 0, Nx(1), 0, Nx(2), 0;
                B.row(1) << 0, Ny(0), 0, Ny(1), 0, Ny(2);
                B.row(2) << Ny(0), Nx(0), Ny(1), Nx(1), Ny(2), Nx(2);

                Nmat = Eigen::MatrixXd::Zero(2, 6);
                Nmat.row(0) << meshElement.shapeFunction2D.N(i,0), 0, meshElement.shapeFunction2D.N(i,1), 0, meshElement.shapeFunction2D.N(i,2), 0;
                Nmat.row(1) << 0, meshElement.shapeFunction2D.N(i,0), 0, meshElement.shapeFunction2D.N(i,1), 0, meshElement.shapeFunction2D.N(i,2);
            }
            // 4 node, linear quadrilateral element {2D}
            else if (meshElement.elemType == 3)
            {
                B = Eigen::MatrixXd::Zero(3, 8);
                B.row(0) << Nx(0), 0, Nx(1), 0, Nx(2), 0, Nx(3), 0;
                B.row(1) << 0, Ny(0), 0, Ny(1), 0, Ny(2), 0, Ny(3);
                B.row(2) << Ny(0), Nx(0), Ny(1), Nx(1), Ny(2), Nx(2), Ny(3), Nx(3);

                Nmat = Eigen::MatrixXd::Zero(2, 8);
                Nmat.row(0) << meshElement.shapeFunction2D.N(i,0), 0, meshElement.shapeFunction2D.N(i,1), 0, meshElement.shapeFunction2D.N(i,2), 0, meshElement.shapeFunction2D.N(i,3), 0;
                Nmat.row(1) << 0, meshElement.shapeFunction2D.N(i,0), 0, meshElement.shapeFunction2D.N(i,1), 0, meshElement.shapeFunction2D.N(i,2), 0, meshElement.shapeFunction2D.N(i,3);

            }
            // 6 node, second order triangle element {2D}
            else if (meshElement.elemType == 9)
            {
                B = Eigen::MatrixXd::Zero(3, 12);
                B.row(0) << Nx(0), 0, Nx(1), 0, Nx(2), 0, Nx(3), 0, Nx(4), 0, Nx(5), 0;
                B.row(1) << 0, Ny(0), 0, Ny(1), 0, Ny(2), 0, Ny(3), 0, Ny(4), 0, Ny(5);
                B.row(2) << Ny(0), Nx(0), Ny(1), Nx(1), Ny(2), Nx(2), Ny(3), Nx(3), Ny(4), Nx(4), Ny(5), Nx(5);

                Nmat = Eigen::MatrixXd::Zero(2, 12);
                Nmat.row(0) << meshElement.shapeFunction2D.N(i,0), 0, meshElement.shapeFunction2D.N(i,1), 0, meshElement.shapeFunction2D.N(i,2), 0, meshElement.shapeFunction2D.N(i,3), 0, meshElement.shapeFunction2D.N(i,4), 0, meshElement.shapeFunction2D.N(i,5), 0;
                Nmat.row(1) << 0, meshElement.shapeFunction2D.N(i,0), 0, meshElement.shapeFunction2D.N(i,1), 0, meshElement.shapeFunction2D.N(i,2), 0, meshElement.shapeFunction2D.N(i,3), 0, meshElement.shapeFunction2D.N(i,4), 0, meshElement.shapeFunction2D.N(i,5);
            }
            // Quadrilateral element, 2nd order {9 noded, quadratic interpolation} (2D)
            else if (meshElement.elemType == 10)
            {
                B = Eigen::MatrixXd::Zero(3, 18);
                B.row(0) << Nx(0), 0, Nx(1), 0, Nx(2), 0, Nx(3), 0, Nx(4), 0, Nx(5), 0, Nx(6), 0, Nx(7), 0, Nx(8), 0 ;
                B.row(1) << 0, Ny(0), 0, Ny(1), 0, Ny(2), 0, Ny(3), 0, Ny(4), 0, Ny(5), 0, Ny(6), 0, Ny(7), 0, Ny(8);
                B.row(2) << Ny(0), Nx(0), Ny(1), Nx(1), Ny(2), Nx(2), Ny(3), Nx(3), Ny(4), Nx(4), Ny(5), Nx(5), Ny(6), Nx(6), Ny(7), Nx(7), Ny(8), Nx(8);

                Nmat = Eigen::MatrixXd::Zero(2, 18);
                Nmat.row(0) << meshElement.shapeFunction2D.N(i,0), 0, meshElement.shapeFunction2D.N(i,1), 0, meshElement.shapeFunction2D.N(i,2), 0, meshElement.shapeFunction2D.N(i,3), 0, meshElement.shapeFunction2D.N(i,4), 0, meshElement.shapeFunction2D.N(i,5), 0, meshElement.shapeFunction2D.N(i,6), 0, meshElement.shapeFunction2D.N(i,7), 0, meshElement.shapeFunction2D.N(i,8), 0;
                Nmat.row(1) << 0, meshElement.shapeFunction2D.N(i,0), 0, meshElement.shapeFunction2D.N(i,1), 0, meshElement.shapeFunction2D.N(i,2), 0, meshElement.shapeFunction2D.N(i,3), 0, meshElement.shapeFunction2D.N(i,4), 0, meshElement.shapeFunction2D.N(i,5), 0, meshElement.shapeFunction2D.N(i,6), 0, meshElement.shapeFunction2D.N(i,7), 0, meshElement.shapeFunction2D.N(i,8);
            }
            // 8 Noded Hexahedral Linear element {3D}
            else if (meshElement.elemType == 5)
            {
                B = Eigen::MatrixXd::Zero(6,24);
                B.row(0) << Nx(0), 0, 0, Nx(1), 0, 0, Nx(2), 0, 0, Nx(3), 0, 0, Nx(4), 0, 0, Nx(5), 0, 0, Nx(6), 0, 0, Nx(7), 0, 0;
                B.row(1) << 0, Ny(0), 0, 0, Ny(1), 0, 0, Ny(2), 0, 0, Ny(3), 0, 0, Ny(4), 0, 0, Ny(5), 0, 0, Ny(6), 0, 0, Ny(7), 0;
                B.row(2) << 0, 0, Nz(0), 0, 0, Nz(1), 0, 0, Nz(2), 0, 0, Nz(3), 0, 0, Nz(4), 0, 0, Nz(5), 0, 0, Nz(6), 0, 0, Nz(7);
                B.row(3) << Ny(0), Nx(0), 0, Ny(1), Nx(1), 0, Ny(2), Nx(2), 0, Ny(3), Nx(3), 0, Ny(4), Nx(4), 0, Ny(5), Nx(5), 0, Ny(6), Nx(6), 0, Ny(7), Nx(7), 0;
                B.row(4) << 0, Nz(0), Ny(0), 0, Nz(1), Ny(1), 0, Nz(2), Ny(2), 0, Nz(3), Ny(3), 0, Nz(4), Ny(4), 0, Nz(5), Ny(5), 0, Nz(6), Ny(6), 0, Nz(7), Ny(7);
                B.row(5) << Nz(0), 0, Nx(0), Nz(1), 0, Nx(1), Nz(2), 0, Nx(2), Nz(3), 0, Nx(3), Nz(4), 0, Nx(4), Nz(5), 0, Nx(5), Nz(6), 0, Nx(6), Nz(7), 0, Nx(7);

                Nmat = Eigen::MatrixXd::Zero(3, 24);
                Nmat.row(0) << meshElement.shapeFunction3D.N(i,0), 0, 0, meshElement.shapeFunction3D.N(i,1), 0, 0, meshElement.shapeFunction3D.N(i,2), 0, 0, meshElement.shapeFunction3D.N(i,3), 0, 0, meshElement.shapeFunction3D.N(i,4), 0, 0, meshElement.shapeFunction3D.N(i,5), 0, 0, meshElement.shapeFunction3D.N(i,6), 0, 0, meshElement.shapeFunction3D.N(i,7), 0, 0;
                Nmat.row(1) << 0, meshElement.shapeFunction3D.N(i,0), 0, 0, meshElement.shapeFunction3D.N(i,1), 0, 0, meshElement.shapeFunction3D.N(i,2), 0, 0, meshElement.shapeFunction3D.N(i,3), 0, 0, meshElement.shapeFunction3D.N(i,4), 0, 0, meshElement.shapeFunction3D.N(i,5), 0, 0, meshElement.shapeFunction3D.N(i,6), 0, 0, meshElement.shapeFunction3D.N(i,7), 0;
                Nmat.row(2) << 0, 0, meshElement.shapeFunction3D.N(i,0), 0, 0, meshElement.shapeFunction3D.N(i,1), 0, 0, meshElement.shapeFunction3D.N(i,2), 0, 0, meshElement.shapeFunction3D.N(i,3), 0, 0, meshElement.shapeFunction3D.N(i,4), 0, 0, meshElement.shapeFunction3D.N(i,5), 0, 0, meshElement.shapeFunction3D.N(i,6), 0, 0, meshElement.shapeFunction3D.N(i,7);
            }
            else{
                std::cerr<<"Element Selection Error : Element Type not available in Solver\n";
                exit(-403);
            }

            // Determinant of Jacobian
            double detJ = Jacobian.determinant();
            detJ = abs(detJ);


            if (solverInp.dimension == 2){
                // Calculate Local element stiffness matrix at each gauss point for 2D element
                K_local += (B.transpose() * Dmat * B) * material.thickness * meshElement.shapeFunction2D.WGP[i] * detJ;
            }
            else if (solverInp.dimension == 3){
                // Calculate Local element stiffness matrix at each gauss point for 3D element
                K_local += (B.transpose() * Dmat * B) * meshElement.shapeFunction3D.WGP[i] * detJ;
            }

            // Vector of body forces {if defined }
            Eigen::VectorXd elementBodyForce;

            // If body force is defined, only then calculate contribution to force vector
            if (bodyForceIndex != -404){

                if (solverInp.dimension == 2){
                    elementBodyForce = Eigen::VectorXd::Zero(2);
                    elementBodyForce << boundary.dirichlet[bodyForceIndex].values[0], boundary.dirichlet[bodyForceIndex].values[1];

                    // Calculate element force vector (local) for 2D elements
                    f_local += material.thickness * Nmat.transpose() * elementBodyForce * meshElement.shapeFunction2D.WGP[i] * detJ;
                }
                else if (solverInp.dimension == 3){
                    elementBodyForce = Eigen::VectorXd::Zero(3);
                    elementBodyForce << boundary.dirichlet[bodyForceIndex].values[0], boundary.dirichlet[bodyForceIndex].values[1],
                                        boundary.dirichlet[bodyForceIndex].values[2];

                    // Calculate element force vector (local) for 3D elements
                    f_local += Nmat.transpose() * elementBodyForce * meshElement.shapeFunction3D.WGP[i] * detJ;
                }
            }
        }
    }
}




// Function to implement Matrix assembly for Linear elasticity problem
void assembleMatrixLE(std::vector<Eigen::Triplet<double>> &tripletList, Eigen::VectorXd &f_global, const std::vector<std::vector<double>> &NODE_COORD, const Element &meshElement, int domainElemID, const SolverInp &solverInp, std::vector<LinearElasticMaterial> material, const BoundaryConditions &boundary, int eqnIndex)
{
    // Initialize vector variable to zeros
    f_global.setZero();

   // Element Connectivity Vector, saves element node connectivities from global connectivity vector
    std::vector<unsigned int> Te(meshElement.ElemConnectivity[0].size());

    // Degree of freedom at each node for the equation we are solving.
    int dof = solverInp.equations[eqnIndex].DOF;


    // std::cout<< " Problem Type "<<solverInp.coordinateSystem << std::endl;

    // if (input.problemType == "AXSYM"){
    //     // For machine part moving with a constant angular velocity        
    //     if (input.bForce[0] == 1){
    //         std::cout<<"Centrifugal forces in rotating parts (if not a rotating part give omega = 0).\nEnter value of constant angular velocity (omega) = ";
    //         std::cin>>input.omega;
    //         std::cout<<"Enter material density (rho) per unit volume = ";
    //         std::cin>>input.RHO;
    //     }
    //     else if (input.bForce[1] == 1)
    //     {
    //         std::cout<<"Enter material density (rho) per unit volume = ";
    //         std::cin>>input.RHO;
    //     }
    // }


    // Loop over elements
    for (int ie = 0; ie < meshElement.numElems; ++ie)
    {
        // Element Connectivity
        Te = meshElement.ElemConnectivity[ie];

        // Vector of element global nodal numbers, as per degree of freedom
        std::vector<unsigned int> Te2dof(Te.size() * dof);

        int ct = 0;
        for (int i = 0, n = Te.size(); i < n; i++)
        {
            // Element node numbers converted to global node nubers (global matrix) as per DOF for each node.
            // Index counter variable (ct)
            for (int j = 0; j < dof; j++)
            {
                Te2dof[ct] = dof * (Te[i] - 1) + j;
                ct++;
            }
        }

        // Element Node Coordinates
        Eigen::MatrixXd Element_NC(Te.size(), solverInp.dimension);
        Element_NC.setZero();

        ct = 0;
        for (unsigned int x : Te)
        {
            if (solverInp.dimension == 3)
            {
                Element_NC.row(ct) << NODE_COORD[x-1][0], NODE_COORD[x-1][1], NODE_COORD[x-1][2];
            }
            else if (solverInp.dimension == 2){
                Element_NC.row(ct) << NODE_COORD[x-1][0], NODE_COORD[x-1][1];
            }
            ct++;
        }


        // Initialize local stiffness element matrix (dense) and force vector
        Eigen::MatrixXd K_local(Element_NC.rows()*dof, Element_NC.rows()*dof);
        Eigen::VectorXd f_local(Element_NC.rows()*dof);

        Eigen::MatrixXd hookesTensor;
        if (solverInp.coordinateSystem == "AXIS"){
            hookesTensor = material[0].HookesLawTensor(-1);
        }
        else{
            hookesTensor = material[0].HookesLawTensor(solverInp.equations[eqnIndex].solverEq);
        }


        // Calculate Element stiffness matrix and Force vector
        ElementMatrixLE(K_local, f_local, Element_NC, hookesTensor, meshElement, solverInp, material[0], boundary);


        // Assemble local stiffness matrix into global matrix
        for (int i = 0, ln = K_local.rows(); i < ln; i++)
        {
            for (int j = 0; j < ln; j++){

                // Adding element stiffness matrix elements to tripletList;
                // (row, column, value at that row and column);
                tripletList.push_back(Eigen::Triplet<double> (Te2dof[i], Te2dof[j], K_local(i,j)));
            }
            f_global(Te2dof[i]) += f_local(i);
        }
    }

    // Converting triplet to sparse matrix (adding to sparse matrix)
    // K_global.setFromTriplets(tripletList.begin(), tripletList.end());
}





// Assmeble global matrix for the heat transfer problem.
void assembleMatrixHT(std::vector<Eigen::Triplet<double>> &tripletStiffness, std::vector<Eigen::Triplet<double>> &tripletMass, Eigen::VectorXd &f_global, const std::vector<std::vector<double>> &NODE_COORD, const Element meshElement[], int domainElemID, const SolverInp &solverInp, const std::vector<MaterialThermal> &material, const BoundaryConditions &boundary)
{

    // Element Connectivity Vector,only element node connectivities from global connectivity vector
    std::vector<unsigned int> Te(meshElement[domainElemID].ElemConnectivity[0].size());


    // Loop over elements of domain
    // std::cout<< "Elemental Loop Starts Here" << std::endl;
    for (int ie = 0; ie < meshElement[domainElemID].numElems; ++ie)
    {
        Te = meshElement[domainElemID].ElemConnectivity[ie];

        // Element Node Coordinates
        Eigen::MatrixXd Element_NC(Te.size(), solverInp.dimension);
        Element_NC.setZero();

        int ct = 0;
        for (unsigned int x : Te)
        {
            if (solverInp.dimension == 3)
            {
                Element_NC.row(ct) << NODE_COORD[x-1][0], NODE_COORD[x-1][1], NODE_COORD[x-1][2];
            }
            else{
                Element_NC.row(ct) << NODE_COORD[x-1][0], NODE_COORD[x-1][1];
            }
            ct++;
        }


        // Initialize local stiffness element matrix (dense) and force vector
        Eigen::MatrixXd K_local(Element_NC.rows(), Element_NC.rows());
        Eigen::MatrixXd C_local(Element_NC.rows(), Element_NC.rows());
        Eigen::VectorXd f_local(Element_NC.rows());


        // ##################### To be changed ####################
        // Currently assuming single equation and single material
        // Calculate Element stiffness matrix and Force vector {cartesian coordinate using isoparametric element}
        int elementTagID = solverInp.equations[0].elemTagId;
        ElementMatrixHT(K_local, C_local, f_local, Element_NC, meshElement[elementTagID], solverInp, material[0], boundary);

        // std::cout<<"Element ("<<ie<<") Klocal = \n"<<K_local<<std::endl;


        // Assemble local stiffness matrix into global matrix
        for (int i = 0, ln = K_local.rows(); i < ln; i++)
        {
            for (int j = 0, lj = K_local.cols(); j < lj; j++)
            {
                // Adding element stiffness matrix elements to tripletList;
                // (row, column, value at that row and column);
                tripletStiffness.push_back(Eigen::Triplet<double> (Te[i]-1, Te[j]-1, K_local(i,j)));

                // Mass Matrix contributions
                if (solverInp.isTransient)
                {
                    tripletMass.push_back(Eigen::Triplet<double> (Te[i]-1, Te[j]-1, C_local(i,j)));
                }
            }

            f_global(Te[i]-1) += f_local(i);
        }
    }
}

