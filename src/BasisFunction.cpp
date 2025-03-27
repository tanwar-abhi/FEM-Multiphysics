

#include "BasisFunction.hpp"



BasisFn::BasisFn(){
    // Default constructor
}

BasisFn::~BasisFn(){
    // Default Destructor
}


// Constructor for creating BasisFun Object
BasisFn::BasisFn(int ngp, int eType)
{
    this->NGP = ngp;
    this->Etype = eType;
}


ShapeFn3D::ShapeFn3D()
{
    // Default Constructor
}

ShapeFn3D::~ShapeFn3D()
{
    // Default Destructor
}




std::vector<Eigen::MatrixXd> BasisFn::ShapeFn_1D(int const ngp)
{
    // Element initialized for 1d Gauss Quadratures
    Quad1D GQelement;
    std::vector<Quad1D> GQ;
    GQ = GQelement.GaussQuad(ngp);

    std::vector<Eigen::MatrixXd> result;

    switch (Etype)
    {
        // 2-node line element
        case 1:
        {
            // Shpae function Matrix of size = ngp x 2
            Eigen::MatrixXd N(ngp, 2);

            // Nxi is derivative of N w.r.t xi;
            Eigen::MatrixXd Nxi(ngp, 2);

            for (int i = 0; i < ngp; i++){

                N.row(i) << (1-GQ[i].zgp1)/2 , (1+GQ[i].zgp1)/2;
                Nxi.row(i) << -0.5 , 0.5;
                this->wgp.push_back(GQ[i].wgp);
            }
            result.push_back(N); result.push_back(Nxi);
            break;
        }
        // 3-node second order line element
        case 8:
        {
            // Shpae function Matrix of size = ngp x 3
            Eigen::MatrixXd N(ngp, 3);

            // Nxi is derivative of N w.r.t xi;
            Eigen::MatrixXd Nxi(ngp, 3);

            // N2xi is derivative of N w.r.t xi twice, i.e. d2N/dxi2;
            Eigen::MatrixXd N2xi(ngp, 3);

            for (int i = 0; i < ngp; i++){

                N.row(i) << (GQ[i].zgp1 - 1)*GQ[i].zgp1/2 , 1 - pow(GQ[i].zgp1,2) , (GQ[i].zgp1 + 1) * GQ[i].zgp1/2;
                Nxi.row(i) << GQ[i].zgp1 - 1/2 , -2 * GQ[i].zgp1 , GQ[i].zgp1 + 1/2;
                N2xi.row(i) << 1.0 , -2.0 , 1.0;
                
                this->wgp.push_back(GQ[i].wgp);
            }
            
            result.push_back(N);
            result.push_back(Nxi); result.push_back(N2xi);
            break;
        }
        // Cubic line element (4-node third order edge)
        // case 26:
        // {
        //     Eigen::MatrixXd N(ngp,4), Nxi(ngp,4), N2xi(ngp,4);

        //     for (int i = 0; i < ngp; i++)
        //     {
        //         N.row(i) << (-9.0/16.0)*(GQ[i].zgp1+1.0/3.0)*(GQ[i].zgp1-1.0/3.0)*(GQ[i].zgp1-1) , (27.0/16.0)*(GQ[i].zgp1+1)*(GQ[i].zgp1-1.0/3.0)*(GQ[i].zgp1-1),
        //                     (-27.0/16.0)*(GQ[i].zgp1+1)*(GQ[i].zgp1+1.0/3.0)*(GQ[i].zgp1-1) , (9.0/16.0)*(GQ[i].zgp1+1)*(GQ[i].zgp1+1.0/3.0)*(GQ[i].zgp1-1.0/3.0);

        //     }

        // }

        default:{
            std::cerr<<"Gauss Points ERROR : Invalid number of gauss points selected.\n";
            exit(-403);
            break;
        }
    }
    return result;
}


std::vector<Eigen::MatrixXd> BasisFn::ShapeFn_2D(const int NGP)
{
    Quad2d GQelement;
    std::vector<Quad2d> GQ = GQelement.GaussQuad(NGP, Etype);

    // Vector containing all values of shape functions
    std::vector<Eigen::MatrixXd> result;

    switch (Etype)
    {
        // 3-node triangle (Linear)
        case 2:{
            // Shpae function Matrix of size = ngp x 3
            Eigen::MatrixXd N(NGP, 3);

            // Nxi is derivative of N w.r.t xi; N2xi is double derivative i.e. d2N/dxi2;
            Eigen::MatrixXd Nxi(NGP, 3), N2xi(NGP, 3);

            // Neta is derivative of N w.r.t. eta; N2eta is double derivative i.e. d2N/deta2;
            Eigen::MatrixXd Neta(NGP, 3), N2eta(NGP, 3);

            // N2xieta is partial detivatieve of N w.r.t xi and eta; N2etaxi vice versa
            Eigen::MatrixXd N2xieta(NGP, 3), N2etaxi(NGP, 3);

            for (int i = 0; i < NGP; i++)
            {
                N.row(i) << 1 - GQ[i].zgp1 - GQ[i].zgp2 , GQ[i].zgp1 , GQ[i].zgp2;
                Nxi.row(i) << -1.0, 1.0, 0;
                Neta.row(i) << -1.0, 0, 1.0;

                this->wgp.push_back(GQ[i].wgp);
            }
            N2xi = Eigen::MatrixXd::Zero(NGP,3);
            N2eta = Eigen::MatrixXd::Zero(NGP,3);
            N2xieta = Eigen::MatrixXd::Zero(NGP,3);
            N2etaxi = Eigen::MatrixXd::Zero(NGP,3);

            result.push_back(N);
            result.push_back(Nxi); result.push_back(N2xi);
            result.push_back(Neta); result.push_back(N2eta);
            result.push_back(N2xieta); result.push_back(N2etaxi);
            break;
        }

        // Quadrilateral element (4-node linear)
        case 3:{
            // Shpae function Matrix of size = ngp x 4
            Eigen::MatrixXd N(NGP, 4);

            // Nxi is derivative of N w.r.t xi; N2xi is double derivative i.e. d2N/dxi2;
            Eigen::MatrixXd Nxi(NGP, 4), N2xi(NGP, 4);

            // Neta is derivative of N w.r.t. eta; N2eta is double derivative i.e. d2N/deta2;
            Eigen::MatrixXd Neta(NGP, 4), N2eta(NGP, 4);

            // N2xieta is partial detivatieve of 'N' w.r.t xi and eta; N2etaxi vice versa
            Eigen::MatrixXd N2xieta(NGP, 4), N2etaxi(NGP, 4);

            for (int i = 0; i < NGP; i++){

                N.row(i) << (1-GQ[i].zgp1)*(1-GQ[i].zgp2)/4 , (1+GQ[i].zgp1)*(1-GQ[i].zgp2)/4 , (1+GQ[i].zgp1)*(1+GQ[i].zgp2)/4 , (1-GQ[i].zgp1)*(1+GQ[i].zgp2)/4;
                Nxi.row(i) << (GQ[i].zgp2-1)/4, (1-GQ[i].zgp2)/4, (1+GQ[i].zgp2)/4, -(1+GQ[i].zgp2)/4;
                Neta.row(i) << (GQ[i].zgp1-1)/4 , -(1+GQ[i].zgp1)/4 , (1+GQ[i].zgp1)/4 , (1-GQ[i].zgp1)/4;
                N2xieta.row(i) << 1.0/4.0 , -1.0/4.0 , 1.0/4.0 , -1.0/4.0;
                N2etaxi.row(i) << 1.0/4.0 , -1.0/4.0 , 1.0/4.0 , -1.0/4.0;

                this->wgp.push_back(GQ[i].wgp);
            }
            N2xi.setZero(); N2eta.setZero();
            
            result.push_back(N);
            result.push_back(Nxi); result.push_back(N2xi);
            result.push_back(Neta); result.push_back(N2eta);
            result.push_back(N2xieta); result.push_back(N2etaxi);
            break;
        }

        // 6-node second order triangle
        case 9:{
            // Shpae function Matrix of size = ngp x 6
            Eigen::MatrixXd N(NGP, 6);

            // Nxi is derivative of N w.r.t xi; N2xi is double derivative i.e. d2N/dxi2;
            Eigen::MatrixXd Nxi(NGP, 6), N2xi(NGP, 6);

            // Neta is derivative of N w.r.t. eta; N2eta is double derivative i.e. d2N/deta2;
            Eigen::MatrixXd Neta(NGP, 6), N2eta(NGP, 6);

            // N2xieta is partial detivatieve of 'N' w.r.t xi and eta; N2etaxi vice versa
            Eigen::MatrixXd N2xieta(NGP, 6), N2etaxi(NGP, 6);

            for (int i = 0; i < NGP; i++)
            {
                double kpa = 1 - GQ[i].zgp1 - GQ[i].zgp2;
                int kpx = -1;
                N.row(i) << GQ[i].zgp1 * (2*GQ[i].zgp1 - 1) , GQ[i].zgp2 * (2*GQ[i].zgp2 - 1) , kpa * (2*kpa - 1) , 4 * GQ[i].zgp1 * GQ[i].zgp2 , 4*GQ[i].zgp2*kpa, 4 * kpa*GQ[i].zgp1 ;
                Nxi.row(i) << 4*GQ[i].zgp1 - 1 , 0 , 4*kpa*kpx-kpx , 4*GQ[i].zgp2 , 4*GQ[i].zgp2*kpx , 4*kpx*GQ[i].zgp1+4*kpa;
                Neta.row(i) << 0 , 4*GQ[i].zgp2-1 , (4*kpa*kpx)-kpx , 4*GQ[i].zgp1 , 4*kpa + 4*GQ[i].zgp2*kpx , 4*kpx*GQ[i].zgp1;
                N2xi.row(i) << 4.0 , 0 , 4.0 , 0 , 0 , -8.0;
                N2eta.row(i) << 0 , 4.0 , 4.0 , 0 , -8.0 , 0;
                N2xieta.row(i) << 0 , 0 , 4.0 , 4.0 , -4.0 , -4.0;
                N2etaxi.row(i) << 0 , 0 , 4.0 , 4.0 , -4.0 , -4.0;

                this->wgp.push_back(GQ[i].wgp);
            }

            result.push_back(N);
            result.push_back(Nxi); result.push_back(N2xi);
            result.push_back(Neta); result.push_back(N2eta);
            result.push_back(N2xieta); result.push_back(N2etaxi);
            break;
        }

        // Quadrilateral element (9-noded, second order quadratic interpolation)
        // (4 nodes associated with the vertices, 4 with the edges and 1 with the face)
        case 10:{
            Eigen::MatrixXd N(NGP, 9);
            Eigen::MatrixXd Nxi(NGP, 9), N2xi(NGP, 9);
            Eigen::MatrixXd Neta(NGP, 9), N2eta(NGP, 9);
            Eigen::MatrixXd N2xieta(NGP, 9), N2etaxi(NGP, 9);

            for (int i = 0; i < NGP; i++)
            {
                N.row(i) << GQ[i].zgp1*(GQ[i].zgp1-1)*GQ[i].zgp2*(GQ[i].zgp2-1)/4.0 , GQ[i].zgp1*(GQ[i].zgp1+1)*GQ[i].zgp2*(GQ[i].zgp2-1)/4.0 ,
                            GQ[i].zgp1*(GQ[i].zgp1+1)*GQ[i].zgp2*(GQ[i].zgp2+1)/4.0 , GQ[i].zgp1*(GQ[i].zgp1-1)*GQ[i].zgp2*(GQ[i].zgp2+1)/4.0 ,
                            (1-pow(GQ[i].zgp1, 2)) * GQ[i].zgp2 * (GQ[i].zgp2-1)/2.0 , GQ[i].zgp1*(GQ[i].zgp1+1)*(1-pow(GQ[i].zgp2,2))/2.0 ,
                            (1-pow(GQ[i].zgp1,2))*GQ[i].zgp2*(GQ[i].zgp2+1)/2.0 , GQ[i].zgp1*(GQ[i].zgp1-1)*(1-pow(GQ[i].zgp2,2))/2.0 ,
                            (1-pow(GQ[i].zgp1,2))*(1-pow(GQ[i].zgp2,2));

                Nxi.row(i) << (GQ[i].zgp1-1.0/2.0)*GQ[i].zgp2*(GQ[i].zgp2-1)/2.0 , (GQ[i].zgp1+1.0/2.0)*GQ[i].zgp2*(GQ[i].zgp2-1)/2.0 ,
                            (GQ[i].zgp1+1.0/2.0)*GQ[i].zgp2*(GQ[i].zgp2+1)/2.0 , (GQ[i].zgp1-1.0/2.0)*GQ[i].zgp2*(GQ[i].zgp2+1)/2.0 ,
                            -GQ[i].zgp1*GQ[i].zgp2*(GQ[i].zgp2-1) , (GQ[i].zgp1+1.0/2.0)*(1-pow(GQ[i].zgp2,2)) ,
                            -GQ[i].zgp1*GQ[i].zgp2*(GQ[i].zgp2+1) , (GQ[i].zgp1-1.0/2.0)*(1-pow(GQ[i].zgp2,2)) ,
                            -2*GQ[i].zgp1*(1- pow(GQ[i].zgp2,2));

                Neta.row(i) << GQ[i].zgp1*(GQ[i].zgp1-1)*(GQ[i].zgp2-1.0/2.0)/2.0 , GQ[i].zgp1*(GQ[i].zgp1+1)*(GQ[i].zgp2-1.0/2.0)/2.0 ,
                            GQ[i].zgp1*(GQ[i].zgp1+1)*(GQ[i].zgp2+1.0/2.0)/2.0 , GQ[i].zgp1*(GQ[i].zgp1-1)*(GQ[i].zgp2+1.0/2.0)/2.0 ,
                            (1-pow(GQ[i].zgp1,2))*(GQ[i].zgp2-1.0/2.0) , GQ[i].zgp1*(GQ[i].zgp1+1)*(-GQ[i].zgp2) ,
                            (1-pow(GQ[i].zgp1,2))*(GQ[i].zgp2+1.0/2.0) , GQ[i].zgp1*(GQ[i].zgp1-1)*(-GQ[i].zgp2) ,
                            (1-pow(GQ[i].zgp1,2))*(-2*GQ[i].zgp2);
                
                N2xi.row(i) << -1.0/2.0*GQ[i].zgp2*(GQ[i].zgp2-1)/2.0 , 1.0/2.0*GQ[i].zgp2*(GQ[i].zgp2-1)/2.0 ,
                            1.0/2.0*GQ[i].zgp2*(GQ[i].zgp2+1)/2.0 , -1.0/2.0*GQ[i].zgp2*(GQ[i].zgp2+1)/2.0 ,
                            -GQ[i].zgp2*(GQ[i].zgp2-1) , 1.0/2.0*(1-pow(GQ[i].zgp2,2)) ,
                            -GQ[i].zgp2*(GQ[i].zgp2+1) , -1.0/2.0*(1-pow(GQ[i].zgp2,2)) ,
                            -2*(1-pow(GQ[i].zgp2,2));

                N2eta.row(i) << -GQ[i].zgp1*(GQ[i].zgp1-1)*1.0/4.0 , -GQ[i].zgp1*(GQ[i].zgp1+1)*1.0/4.0 ,
                            GQ[i].zgp1*(GQ[i].zgp1+1)*1.0/4.0 , GQ[i].zgp1*(GQ[i].zgp1-1)*1.0/4.0 ,
                            -(1-pow(GQ[i].zgp1,2))*1.0/2.0 , -GQ[i].zgp1*(GQ[i].zgp1+1) ,
                            (1-pow(GQ[i].zgp1,2))*1.0/2.0 , -GQ[i].zgp1*(GQ[i].zgp1-1) ,
                            -(1-pow(GQ[i].zgp1,2))*2;

                N2xieta.row(i) << (GQ[i].zgp1-1.0/2.0)*(2*GQ[i].zgp2-1)/2.0 , (GQ[i].zgp1+1.0/2.0)*(2*GQ[i].zgp2-1)/2.0 ,
                            (GQ[i].zgp1+1.0/2.0)*(2*GQ[i].zgp2+1)/2.0 , (GQ[i].zgp1-1.0/2.0)*(2*GQ[i].zgp2+1)/2.0 ,
                            -GQ[i].zgp1*(2*GQ[i].zgp2-1) , (GQ[i].zgp1+1.0/2.0)*(-2*GQ[i].zgp2) ,
                            -GQ[i].zgp1*(2*GQ[i].zgp2+1) , (GQ[i].zgp1-1.0/2.0)*(-2*GQ[i].zgp2) ,
                            4*GQ[i].zgp1*GQ[i].zgp2;

                N2etaxi.row(i) << (2*GQ[i].zgp1-1)*(GQ[i].zgp2-1.0/2.0)/2.0 , (2*GQ[i].zgp1+1)*(GQ[i].zgp2-1.0/2.0)/2.0 ,
                            (2*GQ[i].zgp1+1)*(GQ[i].zgp2+1.0/2.0)/2.0 , (2*GQ[i].zgp1-1)*(GQ[i].zgp2+1.0/2.0)/2.0 ,
                            -2*GQ[i].zgp1*(GQ[i].zgp2-1.0/2.0) , (2*GQ[i].zgp1+1)*(-GQ[i].zgp2) ,
                            -2*GQ[i].zgp1*(GQ[i].zgp2+1.0/2.0) , (2*GQ[i].zgp1-1)*(-GQ[i].zgp2) ,
                            4*GQ[i].zgp1*GQ[i].zgp2;
                
                this->wgp.push_back(GQ[i].wgp);
            }

            result.push_back(N);
            result.push_back(Nxi); result.push_back(N2xi);
            result.push_back(Neta); result.push_back(N2eta);
            result.push_back(N2xieta); result.push_back(N2etaxi);
            break;
        }
        default:{
            std::cerr<<"Element SELECTION ERROR : 2D element selected is not available.\n";
            exit(-403);
            break;
        }
    }
    return result;
}


ShapeFn1D::ShapeFn1D(){
    // Default constructor
}

ShapeFn1D::~ShapeFn1D(){
    // Default Destructor for ShapeFn1D
}


// void ShapeFn1D::getShapeFn(const int gp, int ET)
void ShapeFn1D::getShapeFn(int ET)
{
    // 2noded line, linear line element
    if (ET == 1)
    {
        this->NGP = 2;
    }
    // Qudratic, 3 noded line element
    else if (ET == 8)
    {
        this->NGP = 3;
    }

    // this->NGP = gp;
    
    BasisFn BFn(NGP, ET);
    
    std::vector<Eigen::MatrixXd> results;
    results = BFn.ShapeFn_1D(NGP);

    this->eType = BFn.Etype;
    this->WGP = BFn.wgp;

    N = Eigen::MatrixXd::Zero(results[0].rows(), results[0].cols());
    Nxi = Eigen::MatrixXd::Zero(results[1].rows(), results[1].cols());
    N2xi = Eigen::MatrixXd::Zero(results[0].rows(), results[0].cols());

    N = results[0];
    Nxi = results[1];
    if (results.size() == 3)
    {
        N2xi = results[2];
    }

}




void ShapeFn1D::getShapeFn(int ET, int gp)
{

    this->NGP = gp;

    BasisFn BFn(NGP, ET);

    std::vector<Eigen::MatrixXd> results;
    results = BFn.ShapeFn_1D(NGP);

    this->eType = BFn.Etype;
    this->WGP = BFn.wgp;

    // Initialize to zeros
    N = Eigen::MatrixXd::Zero(results[0].rows(), results[0].cols());
    Nxi = Eigen::MatrixXd::Zero(results[1].rows(), results[1].cols());
    N2xi = Eigen::MatrixXd::Zero(results[0].rows(), results[0].cols());

    N = results[0];
    Nxi = results[1];
    if (results.size() == 3)
    {
        N2xi = results[2];
    }
}




ShapeFn2D::ShapeFn2D(){
    // Default constructor for 2D shape function object
}


ShapeFn2D::~ShapeFn2D(){
    // Default constructor for shape function object
}


// Overload function, used for generating 2D shape functions
// void ShapeFn2D::getShapeFn(int gp, int ET)
void ShapeFn2D::getShapeFn(int ET)
{
    // Setting gauss point for basis function
    // this->NGP = gp;
    
    if (ET == 2)
        this->NGP = 3;
    else if (ET == 9)
        this->NGP = 6;
    else if (ET == 3)
        this->NGP = 4;
    else if (ET == 10){
        this->NGP = 9;
    }

    // Create object of Basis Function for required gauss points and elementType
    BasisFn BFn(NGP, ET);
    
    std::vector<Eigen::MatrixXd> results;
    results = BFn.ShapeFn_2D(NGP);

    // Binds elemnt type and weights for 2D element to object of ShapeFn2D
    this->eType = BFn.Etype;
    this->WGP = BFn.wgp;

    // Bind all shape function and it's derivative to the object of ShapeFn2D
    N = Eigen::MatrixXd::Zero(results[0].rows(), results[0].cols());
    Nxi = Eigen::MatrixXd::Zero(results[1].rows(), results[1].cols());
    N2xi = Eigen::MatrixXd::Zero(results[2].rows(), results[2].cols());
    Neta = Eigen::MatrixXd::Zero(results[3].rows(), results[3].cols());
    N2eta = Eigen::MatrixXd::Zero(results[4].rows(), results[4].cols());
    N2xieta = Eigen::MatrixXd::Zero(results[5].rows(), results[5].cols());
    N2etaxi = Eigen::MatrixXd::Zero(results[6].rows(), results[6].cols());

    N = results[0];
    Nxi = results[1]; N2xi = results[2];
    Neta = results[3]; N2eta = results[4];
    N2xieta = results[5]; N2etaxi = results[6];
}






void ShapeFn2D::getShapeFn(int ET, int gp)
{
    // Setting gauss point for basis function
    this->NGP = gp;

    // Create object of Basis Function for required gauss points and elementType
    BasisFn BFn(NGP, ET);
    
    std::vector<Eigen::MatrixXd> results;
    results = BFn.ShapeFn_2D(NGP);

    // Binds elemnt type and weights for 2D element to object of ShapeFn2D
    this->eType = BFn.Etype;
    this->WGP = BFn.wgp;

    // Bind all shape function and it's derivative to the object of ShapeFn2D
    N = Eigen::MatrixXd::Zero(results[0].rows(), results[0].cols());
    Nxi = Eigen::MatrixXd::Zero(results[1].rows(), results[1].cols());
    N2xi = Eigen::MatrixXd::Zero(results[2].rows(), results[2].cols());
    Neta = Eigen::MatrixXd::Zero(results[3].rows(), results[3].cols());
    N2eta = Eigen::MatrixXd::Zero(results[4].rows(), results[4].cols());
    N2xieta = Eigen::MatrixXd::Zero(results[5].rows(), results[5].cols());
    N2etaxi = Eigen::MatrixXd::Zero(results[6].rows(), results[6].cols());

    N = results[0];
    Nxi = results[1]; N2xi = results[2];
    Neta = results[3]; N2eta = results[4];
    N2xieta = results[5]; N2etaxi = results[6];
}




void ShapeFn3D::getShapeFn(int elemType, int ngp)
{
    // Number of gauss points
    this->NGP = ngp;
    eType = elemType;

    Quad3D gaussQuadElement;
    std::vector<Quad3D> GQ;

    GQ = gaussQuadElement.GaussQuad(ngp, elemType);

    switch (elemType)
    {
        // 8-node hexahedral element (Solid {stiff} element 8 integration points)
        case 5:
        {
            // Initialize all shape functions variables
            N = Eigen::MatrixXd::Zero(ngp, 8); 

            // Nxi,Neta = Differential of N w.r.t xi and eta; 
            // N2xi, N2eta = double differential of shape function w.r.t. xi and eta {natural coordinates}
            Nxi = Eigen::MatrixXd::Zero(ngp, 8); N2xi = Eigen::MatrixXd::Zero(ngp, 8); 
            Neta = Eigen::MatrixXd::Zero(ngp, 8); N2eta = Eigen::MatrixXd::Zero(ngp, 8);
            Nzeta = Eigen::MatrixXd::Zero(ngp, 8); N2zeta = Eigen::MatrixXd::Zero(ngp, 8);

            // Double derivatives of shape function w.r.t xi, eta, zeta
            N2xieta = Eigen::MatrixXd::Zero(ngp, 8); N2XiZeta = Eigen::MatrixXd::Zero(ngp, 8);
            N2etaxi = Eigen::MatrixXd::Zero(ngp, 8); N2EtaZeta = Eigen::MatrixXd::Zero(ngp, 8);
            N2ZetaXi  = Eigen::MatrixXd::Zero(ngp, 8); N2ZetaEta = Eigen::MatrixXd::Zero(ngp, 8);

            for (int i = 0; i < ngp; i++)
            {
                N.row(i) << (1.0/8.0)*(1-GQ[i].zgp1)*(1-GQ[i].zgp2)*(1-GQ[i].zgp3), (1.0/8.0)*(1+GQ[i].zgp1)*(1-GQ[i].zgp2)*(1-GQ[i].zgp3),
                            (1.0/8.0)*(1+GQ[i].zgp1)*(1+GQ[i].zgp2)*(1-GQ[i].zgp3), (1.0/8.0)*(1-GQ[i].zgp1)*(1+GQ[i].zgp2)*(1-GQ[i].zgp3),
                            (1.0/8.0)*(1-GQ[i].zgp1)*(1-GQ[i].zgp2)*(1+GQ[i].zgp3), (1.0/8.0)*(1+GQ[i].zgp1)*(1-GQ[i].zgp2)*(1+GQ[i].zgp3),
                            (1.0/8.0)*(1+GQ[i].zgp1)*(1+GQ[i].zgp2)*(1+GQ[i].zgp3), (1.0/8.0)*(1-GQ[i].zgp1)*(1+GQ[i].zgp2)*(1+GQ[i].zgp3);

                Nxi.row(i) << (-1.0/8.0)*(1-GQ[i].zgp2)*(1-GQ[i].zgp3), (1.0/8.0)*(1-GQ[i].zgp2)*(1-GQ[i].zgp3), (1.0/8.0)*(1+GQ[i].zgp2)*(1-GQ[i].zgp3),
                            (-1.0/8.0)*(1+GQ[i].zgp2)*(1-GQ[i].zgp3), (-1.0/8.0)*(1-GQ[i].zgp2)*(1+GQ[i].zgp3), (1.0/8.0)*(1-GQ[i].zgp2)*(1+GQ[i].zgp3),
                            (1.0/8.0)*(1+GQ[i].zgp2)*(1+GQ[i].zgp3), (-1.0/8.0)*(1+GQ[i].zgp2)*(1+GQ[i].zgp3);

                Neta.row(i) << (-1.0/8.0)*(1-GQ[i].zgp1)*(1-GQ[i].zgp3), (-1.0/8.0)*(1+GQ[i].zgp1)*(1-GQ[i].zgp3), (1.0/8.0)*(1+GQ[i].zgp1)*(1-GQ[i].zgp3),
                                (1.0/8.0)*(1-GQ[i].zgp1)*(1-GQ[i].zgp3), (-1.0/8.0)*(1-GQ[i].zgp1)*(1+GQ[i].zgp3), (-1.0/8.0)*(1+GQ[i].zgp1)*(1+GQ[i].zgp3),
                                (1.0/8.0)*(1+GQ[i].zgp1)*(1+GQ[i].zgp3), (1.0/8.0)*(1-GQ[i].zgp1)*(1+GQ[i].zgp3);

                Nzeta.row(i) << (-1.0/8.0)*(1-GQ[i].zgp1)*(1-GQ[i].zgp2), (-1.0/8.0)*(1+GQ[i].zgp1)*(1-GQ[i].zgp2), (-1.0/8.0)*(1+GQ[i].zgp1)*(1+GQ[i].zgp2),
                                (-1.0/8.0)*(1-GQ[i].zgp1)*(1+GQ[i].zgp2), (1.0/8.0)*(1-GQ[i].zgp1)*(1-GQ[i].zgp2), (1.0/8.0)*(1+GQ[i].zgp1)*(1-GQ[i].zgp2),
                                (1.0/8.0)*(1+GQ[i].zgp1)*(1+GQ[i].zgp2), (1.0/8.0)*(1-GQ[i].zgp1)*(1+GQ[i].zgp2);

                N2xieta.row(i) << (1.0/8.0)*(1-GQ[i].zgp3), (-1.0/8.0)*(1-GQ[i].zgp3), (1.0/8.0)*(1-GQ[i].zgp3), (-1.0/8.0)*(1-GQ[i].zgp3),
                                (1.0/8.0)*(1+GQ[i].zgp3), (-1.0/8.0)*(1+GQ[i].zgp3), (1.0/8.0)*(1+GQ[i].zgp3), (-1.0/8.0)*(1+GQ[i].zgp3);

                N2etaxi.row(i) << N2xieta.row(i);

                N2XiZeta.row(i) << (1.0/8.0)*(1-GQ[i].zgp2), (-1.0/8.0)*(1-GQ[i].zgp2), (-1.0/8.0)*(1+GQ[i].zgp2), (1.0/8.0)*(1+GQ[i].zgp2),
                                    (-1.0/8.0)*(1-GQ[i].zgp2), (1.0/8.0)*(1-GQ[i].zgp2), (1.0/8.0)*(1+GQ[i].zgp2), (-1.0/8.0)*(1+GQ[i].zgp2);

                N2ZetaXi.row(i) << N2XiZeta.row(i);

                N2EtaZeta.row(i) << (1.0/8.0)*(1-GQ[i].zgp1), (1.0/8.0)*(1+GQ[i].zgp1), (-1.0/8.0)*(1+GQ[i].zgp1), (-1.0/8.0)*(1-GQ[i].zgp1),
                                    (-1.0/8.0)*(1-GQ[i].zgp1), (-1.0/8.0)*(1+GQ[i].zgp1), (1.0/8.0)*(1+GQ[i].zgp1), (1.0/8.0)*(1-GQ[i].zgp1);

                N2ZetaEta.row(i) << N2ZetaEta.row(i);

                this->WGP.push_back(GQ[i].wgp);
            }
            break;
        }
        default:{
            std::cerr<<"Gauss Points ERROR : Invalid number of gauss points selected.\n";
            exit(-403);
            break;
        }
    }

}