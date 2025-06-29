
#include "bForce_Traction.hpp"


// Initialize vector of boundary force terms
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




// Imposing neuman boundary forces term in global force vector
void ImposeTraction(ShapeFn1D BasisFn1D, Eigen::VectorXd &fGlobal, std::vector<std::vector<double>> NODE_COORD,std::vector<std::vector<unsigned int>> T, readMesh mesh, Parameters_LE input)
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
        Eigen::MatrixXd Xe(Te.size(), mesh.DIM);
        Xe.setZero();
        int ct = 0;

        for (unsigned int x : Te)
        {
            if (mesh.DIM == 3)
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
            for (int j = 0; j < BasisFn1D.NGP; j++)
            {
                Eigen::MatrixXd Nmat(2, nen*2);
                Nmat.setZero();
                if (nen == 2){
                    Nmat.row(0) << BasisFn1D.N(j,0), 0, BasisFn1D.N(j,1), 0;
                    Nmat.row(1) << 0, BasisFn1D.N(j,0), 0, BasisFn1D.N(j,1);
                }
                else if(nen == 3){
                    Nmat.row(0) << BasisFn1D.N(j,0), 0, BasisFn1D.N(j,1), 0, BasisFn1D.N(j,2), 0;
                    Nmat.row(1) << 0, BasisFn1D.N(j,0), 0, BasisFn1D.N(j,1), 0, BasisFn1D.N(j,2);
                }

                // Boundary element traction force calculation.
                fe_tr += input.thickness * (Nmat.transpose() * (trF/LEdge) * BasisFn1D.WGP[j] * le/2);
                // fe_tr += input.thickness * (Nmat.transpose() * (trF/le) * BasisFn1D.WGP[j] * le/2);
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