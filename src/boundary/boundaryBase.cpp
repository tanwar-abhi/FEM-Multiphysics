
#include "boundaryBase.hpp"


BoundaryBase::BoundaryBase()
{

}

BoundaryBase::~BoundaryBase()
{

}

const std::unordered_map<std::string, BoundaryType> boundaryTypeMap = {
    {"CONSTANT", BoundaryType::Constant},
    {"TEMPORAL", BoundaryType::Temporal},
    {"CONSTANTNODALFILEINPUT", BoundaryType::ConstantNodalFileInput},
    {"TEMPORALNODALFILEINPUT", BoundaryType::TemporalNodalFileInput},
    {"UNKNOWN", BoundaryType::Unknown}
};


BoundaryType BoundaryBase::getBoundaryType()
{
    
}