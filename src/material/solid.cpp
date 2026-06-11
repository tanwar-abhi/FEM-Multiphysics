LinearElasticMaterial::LinearElasticMaterial() {

}

LinearElasticMaterial::~LinearElasticMaterial() {

}

LinearElasticMaterial::LinearElasticMaterial(const float& youngModulus, const float& poissonRatio, const double& angularVelocity=0, const double& density)
{
    E = youngModulus;
    NU = poissonRatio;
    RHO = density;
    if (angularVelocity != 0)
    {
        omega = angularVelocity;
    }
}

// Read material data from JSON input files
std::vector<LinearElasticMaterial> LinearElasticMaterial::readMaterialInputs()
{
    LinearElasticMaterial material;
    std::vector <LinearElasticMaterial> materialVector;

    const std::string solverInputJsonFile = _directoryInput + "solver.json";
    std::ifstream mat_text(solverInputJsonFile);
    Json::Value mat_root;
    Json::Reader mat_reader;
    bool parsingSuccessful = mat_reader.parse(mat_text, mat_root);

    if(!parsingSuccessful)
    {
	    std::cout << "Material Inputs Error : Error parsing the string" << std::endl;
        exit(-404);
    }

    const Json::Value inp_materialProp = mat_root["materialProperty"];

    // Temporary Variable used for reading inputs
    std::string tmp;

    for(int index = 0; index < inp_materialProp.size(); index++)
    {
        tmp = inp_materialProp[index]["type"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

        tmp = inp_materialProp[index]["name"].asString();
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        material.name = tmp;

        material.RHO = inp_materialProp[index]["rho"].asDouble();
        material.NU = inp_materialProp[index]["nu"].asFloat();
        material.thickness = inp_materialProp[index]["thickness"].asFloat();
        material.E = inp_materialProp[index]["youngsMod"].asFloat();
        material.omega = inp_materialProp[index]["omega"].asDouble();

        // Eliminate user input error for thickness
        if (material.thickness <= 0)
        {
            material.thickness = 1;
        }

        materialVector.push_back(material);
    }

    return materialVector;
}

// Calculate constitutive equation Hookes Law Tensor for equation
Eigen::MatrixXd LinearElasticMaterial::HookesLawTensor(int equationType)
{
    // Initializing Linear Elasticity matrix
    Eigen::MatrixXd DMatrix;

    // Plane Stress
    if (equationType == 1)
    {
        DMatrix = Eigen::MatrixXd::Zero(3, 3);
        double term = E/(1 - pow(NU,2));
        DMatrix.row(0) << 1, NU, 0;
        DMatrix.row(1) << NU, 1, 0;
        DMatrix.row(2) << 0, 0, (1-NU)/2;
        DMatrix *= term;
    }
    // Plane Strain
    else if (equationType == 13)
    {
        DMatrix = Eigen::MatrixXd::Zero(3, 3);
        double term = E / ((1+NU)*(1-2*NU));
        DMatrix.row(0) << (1-NU), NU, 0;
        DMatrix.row(1) << NU, (1-NU) , 0;
        DMatrix.row(2) << 0, 0, (1-2*NU)/2;
        DMatrix *= term;
    }
    // Axisymmetric
    else if (equationType == -1)
    {
        DMatrix = Eigen::MatrixXd::Zero(4, 4);
        double term = E / ((1+NU)*(1-2*NU));
        DMatrix.row(0) << 1-NU , NU , NU , 0;
        DMatrix.row(1) << NU , 1-NU , NU , 0;
        DMatrix.row(2) << NU , NU , 1-NU , 0;
        DMatrix.row(3) << 0 , 0 , 0 , (1-2*NU)/2.0;
        DMatrix *= term;
    }
    // 3D Linear elastic
    else if (equationType == 14)
    {
        DMatrix = Eigen::MatrixXd::Zero(6, 6);
        double term = E*(1-NU)/((1+NU)*(1-2*NU));
        DMatrix.row(0) << 1, NU/(1-NU), NU/(1-NU), 0, 0, 0;
        DMatrix.row(1) << NU/(1-NU), 1, NU/(1-NU), 0, 0, 0;
        DMatrix.row(2) << NU/(1-NU), NU/(1-NU), 1, 0, 0, 0;
        DMatrix.row(3) << 0, 0, 0, (1-2*NU)/(2*(1-NU)), 0, 0;
        DMatrix.row(4) << 0, 0, 0, 0, (1-2*NU)/(2*(1-NU)), 0;
        DMatrix.row(5) << 0, 0, 0, 0, 0, (1-2*NU)/(2*(1-NU));
        DMatrix *= term;
    }
    else{
        std::cerr<<"Material Error : Invalid linear ealstic Material properties\n";
        exit(-101);
    }

    return DMatrix;
}