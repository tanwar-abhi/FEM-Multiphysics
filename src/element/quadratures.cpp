
#include <cmath>
#include "Quadratures.hpp"

Quad1D::Quad1D()
{
    // Default Constructor
}

Quad1D::~Quad1D(){
    // Default destructor
}

// Get Gauss Points and weights for 1D element; ngaus -> No. of gauss points
std::vector<Quad1D> Quad1D::GaussQuad(const int NGAUS)
{
    // Initializing vector of quadratures
    std::vector<Quad1D> GPoints(NGAUS);

    switch (NGAUS)
    {
        case 1:
            GPoints[0].zgp1 = 0.0;
            GPoints[0].wgp = 2.0;
        break;

        case 2:
            GPoints[0].zgp1 = -1.0/sqrt(3);
            GPoints[1].zgp1 = 1.0/sqrt(3);

            GPoints[0].wgp = 1.0; GPoints[1].wgp = 1.0;
        break;

        case 3:
            GPoints[0].zgp1 = -0.774596669241483;
            GPoints[1].zgp1 =  0.0;
            GPoints[2].zgp1 =  0.774596669241483;

            GPoints[0].wgp = 5.0/9.0;
            GPoints[1].wgp = 8.0/9.0;
            GPoints[2].wgp = 5.0/9.0;
        break;

        case 4:
            GPoints[0].zgp1 = -0.861136311594953;
            GPoints[1].zgp1 = -0.339981043584856;
            GPoints[2].zgp1 =  0.339981043584856;
            GPoints[3].zgp1 =  0.861136311594953;

            GPoints[0].wgp = 0.347854845137454;
            GPoints[1].wgp = 0.652145154862546;
            GPoints[2].wgp = 0.652145154862546;
            GPoints[3].wgp = 0.347854845137454;
        break;

        case 5:
            GPoints[0].zgp1 = -0.906179845938664;
            GPoints[1].zgp1 = -0.538469310105683;
            GPoints[2].zgp1 =  0.0;
            GPoints[3].zgp1 =  0.538469310105683;
            GPoints[4].zgp1 =  0.906179845938664;

            GPoints[0].wgp = 0.236926885056189;
            GPoints[1].wgp = 0.478628670499366;
            GPoints[2].wgp = 0.568888888888889;
            GPoints[3].wgp = 0.478628670499366;
            GPoints[4].wgp = 0.236926885056189;
        break;
        
        case 6:
            GPoints[0].zgp1 = -0.932469514203152;
            GPoints[1].zgp1 = -0.661209386466265;
            GPoints[2].zgp1 = -0.238619186083197;
            GPoints[3].zgp1 =  0.238619186083197;
            GPoints[4].zgp1 =  0.661209386466265;
            GPoints[5].zgp1 =  0.932469514203152;

            GPoints[0].wgp = 0.171324492379170;
            GPoints[1].wgp = 0.360761573048139;
            GPoints[2].wgp = 0.467913934572691;
            GPoints[3].wgp = 0.467913934572691;
            GPoints[4].wgp = 0.360761573048139;
            GPoints[5].wgp = 0.171324492379170;
        break;

        case 7:
            GPoints[0].zgp1 = -0.9491079123427585245261897;
            GPoints[1].zgp1 = -0.7415311855993944398638648;
            GPoints[2].zgp1 = -0.4058451513773971669066064;
            GPoints[3].zgp1 =  0.0;
            GPoints[4].zgp1 =  0.4058451513773971669066064;
            GPoints[5].zgp1 =  0.7415311855993944398638648;
            GPoints[6].zgp1 =  0.9491079123427585245261897;

            GPoints[0].wgp = 0.1294849661688696932706114;
            GPoints[1].wgp = 0.2797053914892766679014678;
            GPoints[2].wgp = 0.3818300505051189449503698;
            GPoints[3].wgp = 0.4179591836734693877551020;
            GPoints[4].wgp = 0.3818300505051189449503698;
            GPoints[5].wgp = 0.2797053914892766679014678;
            GPoints[6].wgp = 0.1294849661688696932706114;
        break;
        
        case 8:
            GPoints[0].zgp1 = -0.96028986;
            GPoints[1].zgp1 = -0.79666648;
            GPoints[2].zgp1 = -0.52553241;
            GPoints[3].zgp1 = -0.18343464;
            GPoints[4].zgp1 =  0.18343464;
            GPoints[5].zgp1 =  0.52553241;
            GPoints[6].zgp1 =  0.79666648;
            GPoints[7].zgp1 =  0.96028986;

            GPoints[0].wgp = 0.10122854;
            GPoints[1].wgp = 0.22238103;
            GPoints[2].wgp = 0.31370665;
            GPoints[3].wgp = 0.36268378;
            GPoints[4].wgp = 0.36268378;
            GPoints[5].wgp = 0.31370665;
            GPoints[6].wgp = 0.22238103;
            GPoints[7].wgp = 0.10122854;
        break;
        
        case 9:
            GPoints[0].wgp = 0.0812743883615744;
            GPoints[1].wgp = 0.1806481606948574;
            GPoints[2].wgp = 0.2606106964029354;
            GPoints[3].wgp = 0.3123470770400029;
            GPoints[4].wgp = 0.3302393550012598;
            GPoints[5].wgp = 0.3123470770400029;
            GPoints[6].wgp = 0.2606106964029354;
            GPoints[7].wgp = 0.1806481606948574;
            GPoints[8].wgp = 0.0812743883615744;

            GPoints[0].zgp1 = -0.9681602395076261;
            GPoints[1].zgp1 = -0.8360311073266358;
            GPoints[2].zgp1 = -0.6133714327005904;
            GPoints[3].zgp1 = -0.3242534234038089;
            GPoints[4].zgp1 = 0.0000000000000000;
            GPoints[5].zgp1 = 0.3242534234038089;
            GPoints[6].zgp1 = 0.6133714327005904;
            GPoints[7].zgp1 = 0.8360311073266358;
            GPoints[8].zgp1 = 0.9681602395076261;
        break;

        case 10:
            GPoints[0].wgp = 0.0666713443086881;
            GPoints[1].wgp = 0.1494513491505806;
            GPoints[2].wgp = 0.2190863625159820;
            GPoints[3].wgp = 0.2692667193099963;
            GPoints[4].wgp = 0.2955242247147529;
            GPoints[5].wgp = 0.2955242247147529;
            GPoints[6].wgp = 0.2692667193099963;
            GPoints[7].wgp = 0.2190863625159820;
            GPoints[8].wgp = 0.1494513491505806;
            GPoints[9].wgp = 0.0666713443086881;

            GPoints[0].zgp1 = -0.9739065285171717;
            GPoints[1].zgp1 = -0.8650633666889845;
            GPoints[2].zgp1 = -0.6794095682990244;
            GPoints[3].zgp1 = -0.4333953941292472;
            GPoints[4].zgp1 = -0.1488743389816312;
            GPoints[5].zgp1 = 0.1488743389816312;
            GPoints[6].zgp1 = 0.4333953941292472;
            GPoints[7].zgp1 = 0.6794095682990244;
            GPoints[8].zgp1 = 0.8650633666889845;
            GPoints[9].zgp1 = 0.9739065285171717;
        break;

        case 11:
            GPoints[0].wgp = 0.0556685671161737;
            GPoints[1].wgp = 0.1255803694649046;
            GPoints[2].wgp = 0.1862902109277343;
            GPoints[3].wgp = 0.2331937645919905;
            GPoints[4].wgp = 0.2628045445102467;
            GPoints[5].wgp = 0.2729250867779006;
            GPoints[6].wgp = 0.2628045445102467;
            GPoints[7].wgp = 0.2331937645919905;
            GPoints[8].wgp = 0.1862902109277343;
            GPoints[9].wgp = 0.1255803694649046;
            GPoints[10].wgp = 0.0556685671161737;
            
            GPoints[0].zgp1 = -0.9782286581460570;
            GPoints[1].zgp1 = -0.8870625997680953;
            GPoints[2].zgp1 = -0.7301520055740494;
            GPoints[3].zgp1 = -0.5190961292068118;
            GPoints[4].zgp1 = -0.2695431559523450;
            GPoints[5].zgp1 = 0.0000000000000000;
            GPoints[6].zgp1 = 0.2695431559523450;
            GPoints[7].zgp1 = 0.5190961292068118;
            GPoints[8].zgp1 = 0.7301520055740494;
            GPoints[9].zgp1 = 0.8870625997680953;
            GPoints[10].zgp1 = 0.9782286581460570;
        break;

        case 12:
            GPoints[0].wgp = 0.0471753363865118;
            GPoints[1].wgp = 0.1069393259953184;
            GPoints[2].wgp = 0.1600783285433462;
            GPoints[3].wgp = 0.2031674267230659;
            GPoints[4].wgp = 0.2334925365383548;
            GPoints[5].wgp = 0.2491470458134028;
            GPoints[6].wgp = 0.2491470458134028;
            GPoints[7].wgp = 0.2334925365383548;
            GPoints[8].wgp = 0.2031674267230659;
            GPoints[9].wgp = 0.1600783285433462;
            GPoints[10].wgp = 0.1069393259953184;
            GPoints[11].wgp = 0.0471753363865118;

            GPoints[0].zgp1 = -0.9815606342467192;
            GPoints[1].zgp1 = -0.9041172563704749;
            GPoints[2].zgp1 = -0.7699026741943047;
            GPoints[3].zgp1 = -0.5873179542866175;
            GPoints[4].zgp1 = -0.3678314989981802;
            GPoints[5].zgp1 = -0.1252334085114689;
            GPoints[6].zgp1 = 0.1252334085114689;
            GPoints[7].zgp1 = 0.3678314989981802;
            GPoints[8].zgp1 = 0.5873179542866175;
            GPoints[9].zgp1 = 0.7699026741943047;
            GPoints[10].zgp1 = 0.9041172563704749;
            GPoints[11].zgp1 = 0.9815606342467192;
        break;

        case 13:
            GPoints[0].wgp = 0.0404840047653159;
            GPoints[1].wgp = 0.0921214998377285;
            GPoints[2].wgp = 0.1388735102197872;
            GPoints[3].wgp = 0.1781459807619457;
            GPoints[4].wgp = 0.2078160475368885;
            GPoints[5].wgp = 0.2262831802628972;
            GPoints[6].wgp = 0.2325515532308739;
            GPoints[7].wgp = 0.2262831802628972;
            GPoints[8].wgp = 0.2078160475368885;
            GPoints[9].wgp = 0.1781459807619457;
            GPoints[10].wgp = 0.1388735102197872;
            GPoints[11].wgp = 0.0921214998377285;
            GPoints[12].wgp = 0.0404840047653159;

            GPoints[0].zgp1 = -0.9841830547185881;
            GPoints[1].zgp1 = -0.9175983992229779;
            GPoints[2].zgp1 = -0.8015780907333099;
            GPoints[3].zgp1 = -0.6423493394403402;
            GPoints[4].zgp1 = -0.4484927510364469;
            GPoints[5].zgp1 = -0.2304583159551348;
            GPoints[6].zgp1 = 0.0000000000000000;
            GPoints[7].zgp1 = 0.2304583159551348;
            GPoints[8].zgp1 = 0.4484927510364469;
            GPoints[9].zgp1 = 0.6423493394403402;
            GPoints[10].zgp1 = 0.8015780907333099;
            GPoints[11].zgp1 = 0.9175983992229779;
            GPoints[12].zgp1 = 0.9841830547185881;
        break;

        case 14:
            GPoints[0].wgp = 0.0351194603317519;
            GPoints[1].wgp = 0.0801580871597602;
            GPoints[2].wgp = 0.1215185706879032;
            GPoints[3].wgp = 0.1572031671581935;
            GPoints[4].wgp = 0.1855383974779378;
            GPoints[5].wgp = 0.2051984637212956;
            GPoints[6].wgp = 0.2152638534631578;
            GPoints[7].wgp = 0.2152638534631578;
            GPoints[8].wgp = 0.2051984637212956;
            GPoints[9].wgp = 0.1855383974779378;
            GPoints[10].wgp = 0.1572031671581935;
            GPoints[11].wgp = 0.1215185706879032;
            GPoints[12].wgp = 0.0801580871597602;
            GPoints[13].wgp = 0.0351194603317519;

            GPoints[0].zgp1 = -0.9862838086968123;
            GPoints[1].zgp1 = -0.9284348836635735;
            GPoints[2].zgp1 = -0.8272013150697650;
            GPoints[3].zgp1 = -0.6872929048116855;
            GPoints[4].zgp1 = -0.5152486363581541;
            GPoints[5].zgp1 = -0.3191123689278897;
            GPoints[6].zgp1 = -0.1080549487073437;
            GPoints[7].zgp1 = 0.1080549487073437;
            GPoints[8].zgp1 = 0.3191123689278897;
            GPoints[9].zgp1 = 0.5152486363581541;
            GPoints[10].zgp1 = 0.6872929048116855;
            GPoints[11].zgp1 = 0.8272013150697650;
            GPoints[12].zgp1 = 0.9284348836635735;
            GPoints[13].zgp1 = 0.9862838086968123;
        break;

        case 15:
            GPoints[0].wgp = 0.0307532419961173;
            GPoints[1].wgp = 0.0703660474881081;
            GPoints[2].wgp = 0.1071592204671719;
            GPoints[3].wgp = 0.1395706779261543;
            GPoints[4].wgp = 0.1662692058169939;
            GPoints[5].wgp = 0.1861610000155622;
            GPoints[6].wgp = 0.1984314853271116;
            GPoints[7].wgp = 0.2025782419255613;
            GPoints[8].wgp = 0.1984314853271116;
            GPoints[9].wgp = 0.1861610000155622;
            GPoints[10].wgp = 0.1662692058169939;
            GPoints[11].wgp = 0.1395706779261543;
            GPoints[12].wgp = 0.1071592204671719;
            GPoints[13].wgp = 0.0703660474881081;
            GPoints[14].wgp = 0.0307532419961173;

            GPoints[0].zgp1 = -0.9879925180204854;
            GPoints[1].zgp1 = -0.9372733924007060;
            GPoints[2].zgp1 = -0.8482065834104272;
            GPoints[3].zgp1 = -0.7244177313601701;
            GPoints[4].zgp1 = -0.5709721726085388;
            GPoints[5].zgp1 = -0.3941513470775634;
            GPoints[6].zgp1 = -0.2011940939974345;
            GPoints[7].zgp1 =  0.0000000000000000;
            GPoints[8].zgp1 =  0.2011940939974345;
            GPoints[9].zgp1 =  0.3941513470775634;
            GPoints[10].zgp1 = 0.5709721726085388;
            GPoints[11].zgp1 = 0.7244177313601701;
            GPoints[12].zgp1 = 0.8482065834104272;
            GPoints[13].zgp1 = 0.9372733924007060;
            GPoints[14].zgp1 = 0.9879925180204854;            
        break;
        
        default:{
            std::cerr<<"Gauss Points SELECTION ERROR : Invalid number of gauss points selected (1D element).\n";
            exit(-403);
            // break;
        }
    }
    
    return GPoints;
}


Quad2d::Quad2d(){
    // Default Constructor
}


Quad2d::~Quad2d(){
    // Default Destructor
}

// Gauss Points for 2D 
std::vector<Quad2d> Quad2d::GaussQuad(const int NGP, int eType)
{
    std::vector<Quad2d> GPoints2D(NGP);
    // Quadrilateral element
    if (eType == 3 || eType == 10)
    {
        int ng = ceil(sqrt(NGP));

        // Element to get 1D gauss quadratures used to calculate Quadrilateral Gauss points.
        Quad1D GQ1D;
        std::vector<Quad1D> GPVector = GQ1D.GaussQuad(ng);
        // Counters for index of gauss points
        int ct = 0;

        for (int i = 0; i < ng; i++)
        {
            for (int j = 0; j < ng; j++)
            {
                GPoints2D[ct].zgp2 = GPVector[i].zgp1;
                GPoints2D[ct].zgp1 = GPVector[j].zgp1;
                GPoints2D[ct].wgp = GPVector[j].wgp * GPVector[i].wgp;
                ct++;
            }
        }
    }

    // Gauss Points for Triangle element
    else{
        switch (NGP)
        {
            case 1:
            {
                GPoints2D[0].zgp1 = 1.0/3.0;
                GPoints2D[0].zgp2 = 1.0/3.0;
                GPoints2D[0].wgp = 0.5 * 1.0/3.0;
                break;
            }
            
            case 3:
            {
                GPoints2D[0].zgp1 = 1.0/6.0; GPoints2D[0].zgp2 = 1.0/6.0;
                GPoints2D[1].zgp1 = 4.0/6.0; GPoints2D[1].zgp2 = 1.0/6.0;
                GPoints2D[2].zgp1 = 1.0/6.0; GPoints2D[2].zgp2 = 4.0/6.0;

                GPoints2D[0].wgp = 0.5*1.0/3.0;
                GPoints2D[1].wgp = 0.5*1.0/3.0;
                GPoints2D[2].wgp = 0.5*1.0/3.0;
                break;
            }

            case 4:
            {
                GPoints2D[0].zgp1 = 1.0/3.0; GPoints2D[0].zgp2 = 1.0/3.0;
                GPoints2D[1].zgp1 = 0.6; GPoints2D[1].zgp2 = 0.2;
                GPoints2D[2].zgp1 = 0.2; GPoints2D[2].zgp2 = 0.6;
                GPoints2D[3].zgp1 = 0.2; GPoints2D[3].zgp2 = 0.2;

                GPoints2D[0].wgp = 0.5*(-27.0/48.0);
                GPoints2D[1].wgp = 0.5*25.0/48.0;
                GPoints2D[2].wgp = 0.5*25.0/48.0;
                GPoints2D[3].wgp = 0.5*25.0/48.0;
                break;
            }

            case 6:
            {
                GPoints2D[0].zgp1 = 0.10810301816807022736; GPoints2D[0].zgp2 = 0.44594849091596488632;
                GPoints2D[1].zgp1 = 0.44594849091596488632; GPoints2D[1].zgp2 = 0.10810301816807022736;
                GPoints2D[2].zgp1 = 0.44594849091596488632; GPoints2D[2].zgp2 = 0.44594849091596488632;
                GPoints2D[3].zgp1 = 0.81684757298045851308; GPoints2D[3].zgp2 = 0.09157621350977074346;
                GPoints2D[4].zgp1 = 0.09157621350977074346; GPoints2D[4].zgp2 = 0.81684757298045851308;
                GPoints2D[5].zgp1 = 0.09157621350977074346; GPoints2D[5].zgp2 = 0.09157621350977074346;

                GPoints2D[0].wgp = 0.5 * 0.22338158967801146570;
                GPoints2D[1].wgp = 0.5 * 0.22338158967801146570;
                GPoints2D[2].wgp = 0.5 * 0.22338158967801146570;
                GPoints2D[3].wgp = 0.5 * 0.10995174365532186764;
                GPoints2D[4].wgp = 0.5 * 0.10995174365532186764;
                GPoints2D[5].wgp = 0.5 * 0.10995174365532186764;
                break;
            }
            
            case 7:
            {
                GPoints2D[0].zgp1 = 1.0/3.0;                GPoints2D[0].zgp2 = 1.0/3.0;
                GPoints2D[1].zgp1 = 0.79742698535308732240; GPoints2D[1].zgp2 = 0.10128650732345633880;
                GPoints2D[2].zgp1 = 0.10128650732345633880; GPoints2D[2].zgp2 = 0.79742698535308732240;
                GPoints2D[3].zgp1 = 0.10128650732345633880; GPoints2D[3].zgp2 = 0.10128650732345633880;
                GPoints2D[4].zgp1 = 0.05971587178976982045; GPoints2D[4].zgp2 = 0.47014206410511508977;
                GPoints2D[5].zgp1 = 0.47014206410511508977; GPoints2D[5].zgp2 = 0.05971587178976982045;
                GPoints2D[6].zgp1 = 0.47014206410511508977; GPoints2D[6].zgp2 = 0.47014206410511508977;

                GPoints2D[0].wgp = 0.5*0.225;
                GPoints2D[1].wgp = 0.5*0.12593918054482715260;
                GPoints2D[2].wgp = 0.5*0.12593918054482715260;
                GPoints2D[3].wgp = 0.5*0.12593918054482715260;
                GPoints2D[4].wgp = 0.5*0.13239415278850618074;
                GPoints2D[5].wgp = 0.5*0.13239415278850618074;
                GPoints2D[6].wgp = 0.5*0.13239415278850618074;
                break;
            }
            
            case 12:
            {
                GPoints2D[0].zgp1 = 0.87382197101699554332;    GPoints2D[0].zgp2 = 0.06308901449150222834;
                GPoints2D[1].zgp1 = 0.06308901449150222834;    GPoints2D[1].zgp2 = 0.87382197101699554332;
                GPoints2D[2].zgp1 = 0.06308901449150222834;    GPoints2D[2].zgp2 = 0.06308901449150222834;
                GPoints2D[3].zgp1 = 0.50142650965817915742;    GPoints2D[3].zgp2 = 0.24928674517091042129;
                GPoints2D[4].zgp1 = 0.24928674517091042129;    GPoints2D[4].zgp2 = 0.50142650965817915742;
                GPoints2D[5].zgp1 = 0.24928674517091042129;    GPoints2D[5].zgp2 = 0.24928674517091042129;
                GPoints2D[6].zgp1 = 0.05314504984481694735;    GPoints2D[6].zgp2 = 0.31035245103378440542;
                GPoints2D[7].zgp1 = 0.31035245103378440542;    GPoints2D[7].zgp2 = 0.05314504984481694735;
                GPoints2D[8].zgp1 = 0.05314504984481694735;    GPoints2D[8].zgp2 = 0.63650249912139864723;
                GPoints2D[9].zgp1 = 0.31035245103378440542;    GPoints2D[9].zgp2 = 0.63650249912139864723;
                GPoints2D[10].zgp1 = 0.63650249912139864723;   GPoints2D[10].zgp2 = 0.05314504984481694735;
                GPoints2D[11].zgp1 = 0.63650249912139864723;   GPoints2D[11].zgp2 = 0.31035245103378440542;

                GPoints2D[0].wgp = 0.5*0.050844906370206816921;
                GPoints2D[1].wgp = 0.5*0.050844906370206816921;
                GPoints2D[2].wgp = 0.5*0.050844906370206816921;
                GPoints2D[3].wgp = 0.5*0.116786275726379366030;
                GPoints2D[4].wgp = 0.5*0.116786275726379366030;
                GPoints2D[5].wgp = 0.5*0.116786275726379366030;
                GPoints2D[6].wgp = 0.5*0.082851075618373575194;
                GPoints2D[7].wgp = 0.5*0.082851075618373575194;
                GPoints2D[8].wgp = 0.5*0.082851075618373575194;
                GPoints2D[9].wgp = 0.5*0.082851075618373575194;
                GPoints2D[10].wgp = 0.5*0.082851075618373575194;
                GPoints2D[11].wgp = 0.5*0.082851075618373575194;
                break;
            }
          
            case 13:
            {
                GPoints2D[0].zgp1 = 0.33333333333333;  GPoints2D[0].zgp2 = 0.33333333333333;
                GPoints2D[1].zgp1 = 0.26034596607904;  GPoints2D[1].zgp2 = 0.26034596607904;
                GPoints2D[2].zgp1 = 0.26034596607904;  GPoints2D[2].zgp2 = 0.47930806784192;
                GPoints2D[3].zgp1 = 0.47930806784192;  GPoints2D[3].zgp2 = 0.26034596607904;
                GPoints2D[4].zgp1 = 0.06513010290222;  GPoints2D[4].zgp2 = 0.06513010290222;
                GPoints2D[5].zgp1 = 0.06513010290222;  GPoints2D[5].zgp2 = 0.86973979419557;
                GPoints2D[6].zgp1 = 0.86973979419557;  GPoints2D[6].zgp2 = 0.06513010290222;
                GPoints2D[7].zgp1 = 0.31286549600487;  GPoints2D[7].zgp2 = 0.63844418856981;
                GPoints2D[8].zgp1 = 0.63844418856981;  GPoints2D[8].zgp2 = 0.04869031542532;
                GPoints2D[9].zgp1 = 0.04869031542532;  GPoints2D[9].zgp2 = 0.31286549600487;
                GPoints2D[10].zgp1 = 0.63844418856981; GPoints2D[10].zgp2 = 0.31286549600487;
                GPoints2D[11].zgp1 = 0.31286549600487; GPoints2D[11].zgp2 = 0.04869031542532;
                GPoints2D[12].zgp1 = 0.04869031542532; GPoints2D[12].zgp2 = 0.63844418856981;

                GPoints2D[0].wgp = 0.5*-0.14957004446768;
                GPoints2D[1].wgp = 0.5* 0.17561525743321;
                GPoints2D[2].wgp = 0.5* 0.17561525743321;
                GPoints2D[3].wgp = 0.5* 0.17561525743321;
                GPoints2D[4].wgp = 0.5* 0.05334723560884;
                GPoints2D[5].wgp = 0.5* 0.05334723560884;
                GPoints2D[6].wgp = 0.5* 0.05334723560884;
                GPoints2D[7].wgp = 0.5* 0.07711376089026;
                GPoints2D[8].wgp = 0.5* 0.07711376089026;
                GPoints2D[9].wgp = 0.5* 0.07711376089026;
                GPoints2D[10].wgp = 0.5* 0.07711376089026;
                GPoints2D[11].wgp = 0.5* 0.07711376089026;
                GPoints2D[12].wgp = 0.5* 0.07711376089026;
                break;
            }
            
            case 16:
            {
                GPoints2D[0].zgp1 = 0.33333333333333;  GPoints2D[0].zgp2 = 0.33333333333333;
                GPoints2D[1].zgp1 = 0.45929258829272;  GPoints2D[1].zgp2 = 0.45929258829272;
                GPoints2D[2].zgp1 = 0.45929258829272;  GPoints2D[2].zgp2 = 0.08141482341455;
                GPoints2D[3].zgp1 = 0.08141482341455;  GPoints2D[3].zgp2 = 0.45929258829272;
                GPoints2D[4].zgp1 = 0.17056930775176;  GPoints2D[4].zgp2 = 0.17056930775176;
                GPoints2D[5].zgp1 = 0.17056930775176;  GPoints2D[5].zgp2 = 0.65886138449648;
                GPoints2D[6].zgp1 = 0.65886138449648;  GPoints2D[6].zgp2 = 0.17056930775176;
                GPoints2D[7].zgp1 = 0.05054722831703;  GPoints2D[7].zgp2 = 0.05054722831703;
                GPoints2D[8].zgp1 = 0.05054722831703;  GPoints2D[8].zgp2 = 0.89890554336594;
                GPoints2D[9].zgp1 = 0.89890554336594;  GPoints2D[9].zgp2 = 0.05054722831703;
                GPoints2D[10].zgp1 = 0.26311282963464; GPoints2D[10].zgp2 = 0.72849239295540;
                GPoints2D[11].zgp1 = 0.72849239295540; GPoints2D[11].zgp2 = 0.00839477740996;
                GPoints2D[12].zgp1 = 0.00839477740996; GPoints2D[12].zgp2 = 0.26311282963464;
                GPoints2D[13].zgp1 = 0.72849239295540; GPoints2D[13].zgp2 = 0.26311282963464;
                GPoints2D[14].zgp1 = 0.26311282963464; GPoints2D[14].zgp2 = 0.00839477740996;
                GPoints2D[15].zgp1 = 0.00839477740996; GPoints2D[15].zgp2 = 0.72849239295540;

                GPoints2D[0].wgp = 0.5 * 0.14431560767779;
                GPoints2D[1].wgp = 0.5 * 0.09509163426728;
                GPoints2D[2].wgp = 0.5 * 0.09509163426728;
                GPoints2D[3].wgp = 0.5 * 0.09509163426728;
                GPoints2D[4].wgp = 0.5 * 0.10321737053472;
                GPoints2D[5].wgp = 0.5 * 0.10321737053472;
                GPoints2D[6].wgp = 0.5 * 0.10321737053472;
                GPoints2D[7].wgp = 0.5 * 0.03245849762320;
                GPoints2D[8].wgp = 0.5 * 0.03245849762320;
                GPoints2D[9].wgp = 0.5 * 0.03245849762320;
                GPoints2D[10].wgp = 0.5 * 0.02723031417443;
                GPoints2D[11].wgp = 0.5 * 0.02723031417443;
                GPoints2D[12].wgp = 0.5 * 0.02723031417443;
                GPoints2D[13].wgp = 0.5 * 0.02723031417443;
                GPoints2D[14].wgp = 0.5 * 0.02723031417443;
                GPoints2D[15].wgp = 0.5 * 0.02723031417443;
                break;
            }
            
            default:{
                std::cerr << "Gauss Points SELECTION ERROR : Invalid number of gauss points selected (2D element).\n";
                exit(-403);
                // break;
            }
        }
    }

    return GPoints2D;
}



Quad3D::Quad3D(){
    // Default constructor for 3D gauss quadratures
}

Quad3D::~Quad3D(){
    // Default destructor of 3D gauss quadratures
}

std::vector<Quad3D> Quad3D::GaussQuad(const int NGP, int elementType)
{
    std::vector<Quad3D> gaussPoints(NGP);

    // Hexahedral element {Linear 8node} :: QUadratic TODO
    if (elementType == 5)
    {
        switch (NGP)
        {
            case 1:
            {
                gaussPoints[0].zgp1 = 0; gaussPoints[0].zgp2 = 0; gaussPoints[0].zgp3 = 0;
                gaussPoints[0].wgp = 8.0;
                break;
            }

            case 6:
            {
                gaussPoints[0].zgp1 = -1.0; gaussPoints[0].zgp2 = 0; gaussPoints[0].zgp3 = 0;
                gaussPoints[1].zgp1 = 1.0; gaussPoints[1].zgp2 = 0; gaussPoints[1].zgp3 = 0;
                gaussPoints[2].zgp1 = 0; gaussPoints[2].zgp2 = -1.0; gaussPoints[2].zgp3 = 0;
                gaussPoints[3].zgp1 = 0; gaussPoints[3].zgp2 = 1.0; gaussPoints[3].zgp3 = 0;
                gaussPoints[4].zgp1 = 0; gaussPoints[4].zgp2 = 0; gaussPoints[4].zgp3 = -1.0;
                gaussPoints[5].zgp1 = 0; gaussPoints[5].zgp2 = 0; gaussPoints[5].zgp3 = 1.0;

                gaussPoints[0].wgp = 4.0/3.0; gaussPoints[1].wgp = 4.0/3.0; gaussPoints[2].wgp = 4.0/3.0;
                gaussPoints[3].wgp = 4.0/3.0; gaussPoints[4].wgp = 4.0/3.0; gaussPoints[5].wgp = 4.0/3.0;
                break;
            }

            case 8:
            {
                gaussPoints[0].zgp1 = -0.5773502691896258; gaussPoints[0].zgp2 = -0.5773502691896258; gaussPoints[0].zgp3 = -0.5773502691896258;
                gaussPoints[1].zgp1 = 0.5773502691896258; gaussPoints[1].zgp2 = -0.5773502691896258; gaussPoints[1].zgp3 = -0.5773502691896258;
                gaussPoints[2].zgp1 = 0.5773502691896258; gaussPoints[2].zgp2 = 0.5773502691896258; gaussPoints[2].zgp3 = -0.5773502691896258;
                gaussPoints[3].zgp1 = -0.5773502691896258; gaussPoints[3].zgp2 = 0.5773502691896258; gaussPoints[3].zgp3 = -0.5773502691896258;
                gaussPoints[4].zgp1 = -0.5773502691896258; gaussPoints[4].zgp2 = -0.5773502691896258; gaussPoints[4].zgp3 = 0.5773502691896258;
                gaussPoints[5].zgp1 = 0.5773502691896258; gaussPoints[5].zgp2 = -0.5773502691896258; gaussPoints[5].zgp3 = 0.5773502691896258;
                gaussPoints[6].zgp1 = 0.5773502691896258; gaussPoints[6].zgp2 = 0.5773502691896258; gaussPoints[6].zgp3 = 0.5773502691896258;
                gaussPoints[7].zgp1 = -0.5773502691896258; gaussPoints[7].zgp2 = 0.5773502691896258; gaussPoints[7].zgp3 = 0.5773502691896258;

                gaussPoints[0].wgp = 1.0; gaussPoints[1].wgp = 1.0; gaussPoints[2].wgp = 1.0; gaussPoints[3].wgp = 1.0;
                gaussPoints[4].wgp = 1.0; gaussPoints[5].wgp = 1.0; gaussPoints[6].wgp = 1.0; gaussPoints[7].wgp = 1.0;
                break;
            }

            default:
            {
                std::cerr<<"ERORR INVALID GAUSS POINT : Invalid number of gauss points selected (3D element).\n";
                exit(-403);
                break;
            }
        }
    }
    return gaussPoints;
}


