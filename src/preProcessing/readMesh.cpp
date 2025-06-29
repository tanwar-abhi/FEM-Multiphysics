
#include "readMesh.hpp"

readMesh::readMesh(){
    // Default constructor
}

readMesh::~readMesh(){
    // Default Destructor
}


// Read Mesh constructor to initialize filename and dimension of each node
readMesh::readMesh(const PartsInp &partsInpObj, int dimension)
{
    // initialize the value of number of nodes to 0
    this->NNodes = 0;
    //    this->partsInp = partsInpObj;
    this->NODE_DIM = dimension;

    // Deep Copy of partsInpObj. 
    // this->partsInpObj = partsInpObj;
    this->partsInpObj.numParts = partsInpObj.numParts;
    this->partsInpObj.part = new Part[partsInpObj.numParts];
    for (int i = 0; i < partsInpObj.numParts; i++){
        this->partsInpObj.part[i] = partsInpObj.part[i];
    }



    // Read mesh nodes
    // std::cout << "started reading mesh nodes" <<  std::endl;
    this->readNodes();
    // std::cout << "completed reading mesh nodes" << std::endl;

    // Implementation of physical field inside each mesh file
    // std::cout<<"Started reading physical element tags\n";
    this->readElementTagsNew();
    // std::cout<<"Done reading physical element tags\n";

}

// Read the elemental tags from the physical field block of the msh file
void readMesh::readElementTagsNew()
{
    // Open and read file.
    for (int i=0; i < this->partsInpObj.numParts; i++)
    {
        std::ifstream meshFile(partsInpObj.part[i].meshFileName);

	    // Flags to initiate begin reading total numbers of nodes and elements and element tags from msh file
	    bool readETag = false;

	    // Read each line
	    std::string Line;

	    // physicalName Tag Id
	    int elemTagId;

	    // physicalName Tag name
	    std::string elemTagName;

	    /*
	     * Open the mesh files
	     */

	    // If File could't be opened
	    if (!meshFile)
	    {
		    std::cerr<<"Error: Unable to read elemement Tags from file\n";
		    std::cout<<"Invalid file or Path/File Name incorrect\n";
		    exit(1);
	    }

	    /*
	     * Read the mesh line by line and check if we have PhysicalName tag
	     * Then read physicalNames tags and their associated ids
	     */
            // Flag to identify line containing total number of element tags
            bool readNElemTags;
	    while(getline(meshFile, Line))
	    {
            // Remove Carriage return if present in the string
            if (!Line.empty() && Line[Line.size() - 1] == '\r')
			    Line.erase(Line.size() - 1);

          if (Line == "$EndPhysicalNames"){
		        meshFile.close();
		        return;
	        }

	        if (Line == "$PhysicalNames"){
		        readETag = true;
                readNElemTags = true;
                continue;
	        }

            if (readETag){
                // Read total number of elements tags (i.e. physical tags) defined in mesh
                if (readNElemTags)
                {
                    NElementalTags = stoi(Line);
                    readNElemTags = false;
                }
                // If length of line < 3 then, it does't contain physical names for nodes
                if (countWords(Line) > 2 )
		        {
			        this->split(Line,elemTagName,elemTagId);
                    std::transform(elemTagName.begin(), elemTagName.end(), elemTagName.begin(), ::toupper);
			        this->elementalTagIds.push_back(elemTagId) ;
			        this->elementalTags.push_back(elemTagName) ;
                    this->elementalTagsMeshName.push_back(partsInpObj.part[i].meshFileName) ;
                    this->elementalTagsPartsId.push_back(i);
                }
	        }
        }
    }
}

// Split string values from the msh file into physical field names and id
void readMesh::split(std::string Line, std::string &elemTagName, int &elemTagId)
{
    std::stringstream ss(Line), ss1;
    std::string value1, tmp;

    ss >> value1 >> elemTagId;
    ss >> elemTagName;
    
    if (elemTagName[0] == '\"')
    {
	    elemTagName = elemTagName.substr(1,elemTagName.length()-1);
	    if (elemTagName[elemTagName.length()-1] == '\"')
	       elemTagName = elemTagName.substr(0,elemTagName.length()-1);

	    ss1 << elemTagName;

	    while(ss>>tmp)
        {
            if (tmp[tmp.length()-1] == '\"')
               tmp = tmp.substr(0,tmp.length()-1);

            ss1<<' '<<tmp;
	    }
    }
    elemTagName = ss1.str();
}


// Split string values from msh file into nodal coordinates vector
void readMesh::pushNodeCoordString(std::string Line, std::vector<std::vector<double>> &NodeCoords, const int DIM)
{
    // Creating string stream of line to split values by using empty space
    std::stringstream ss(Line);
    std::string value1, value2, value3, value4;
    ss >> value1 >> value2 >> value3 >> value4;

    // unsigned int pos = stoul(value1);
    double X = stod(value2), Y = stod(value3), Z = stod(value4);

    // Adding X and Y Nodal coordinates
    std::vector<double> coord;
    coord.push_back(X);
    coord.push_back(Y);
    coord.push_back(Z);

    NodeCoords.push_back(coord);
}




// Read node coordinates and return a 2D vector of nodal coordinates
void readMesh::readNodes()
{
    // Open and read file.
    for (int i=0; i < this->partsInpObj.numParts; i++)
    {
      startingNodeIdMesh.push_back(NNodes);
      std::ifstream meshFile(partsInpObj.part[i].meshFileName);

      // If File could't be opened
      if (!meshFile)
      {
          std::cerr<<"Error: Unable to read Nodes from file\n";
          std::cout<<"Invalid file or Path/File Name incorrect\n";
          exit(1);
      }

      // For Read each line
      std::string eachLine;

      // Read node flag, tell the function to start assigning values to vector elements
      bool readFlag = false;

      while(getline(meshFile, eachLine))
      {
	    if (!eachLine.empty() && eachLine[eachLine.size() - 1] == '\r')
            eachLine.erase(eachLine.size() - 1);

          // "\r" is the carriage return i.e. moves cursor to start of next line
          if (eachLine == "$Nodes" || eachLine == "$ParametricNodes")
          {
              readFlag = true;
              continue;
          }
          else if (eachLine == "$EndNodes" || eachLine == "$EndParametricNodes"){
              readFlag = false;
              break;
          }
        
          if (readFlag)
          {
              if (eachLine.length() > 5)
              {
                  pushNodeCoordString(eachLine, this->Node_Coord, this->NODE_DIM);
              }
          }
      }
    
      meshFile.close();
      this->NNodes = this->NNodes + this->Node_Coord.size();
    }
}

int countWords(std::string str)
{

   int wordCount;

   std::stringstream ss(str);
   std::string word;

   // Check if the string is null
   // or empty then return zero
   if (str.length()==0)
        return 0;

   wordCount = 0;
   while (ss >> word){
       wordCount++;
   }

    return wordCount;
}
