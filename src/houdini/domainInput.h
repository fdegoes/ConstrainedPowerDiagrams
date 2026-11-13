#ifndef PP_HOUDINI_DOMAININPUT_H
#define PP_HOUDINI_DOMAININPUT_H

#include <GA/GA_Types.h>

#include <string>
#include <vector>

class GU_Detail;
class PRM_ChoiceList;
class PRM_Name;

namespace PP
{

class Domain;

// "mode" menu shared by the SOPs that take a clipping domain.
enum DomainMode { MODE_2D = 0, MODE_3D = 1 };

extern PRM_Name       modeParmName;
extern PRM_ChoiceList modeParmMenu;

// Sites of input 0, projected on the working space (the chosen plane in
// 2D, R^3 in 3D).
struct SiteInput
{
    int                    _dim = 0;
    std::vector<GA_Offset> _offsets;  // point offsets, in order
    std::vector<double>    _sites;    // effective sites P + u, dim per site
    std::vector<double>    _weights;  // effective weights w + |u|^2
    std::vector<double>    _shiftSq;  // |u|^2 (in-plane part in 2D)
    std::string            _warning;  // missing attributes
};

// Reads the points of gdp. weightAttr and shiftAttr may be empty
// or missing (then 0). In 2D, only the in-plane part of the shift is used.
void ReadSites(
    const GU_Detail* gdp,
    bool is2D,
    const int* axes,
    const char* weightAttr,
    const char* shiftAttr,
    SiteInput& out);

// Reads the clipping domain of gdp: its triangles (2D, projected on the
// plane) or tetrahedra (3D).
//
// pts         dim doubles per point of gdp
// elems       3 (2D) or 4 (3D) point indices per element
// planeOffset 2D: average height of the triangles along the normal
bool ReadDomain(
    const GU_Detail* gdp,
    bool is2D,
    const int* axes,
    std::vector<double>& pts,
    std::vector<int>& elems,
    double& planeOffset,
    std::string& err);

// Reads the clipping domain of gdp and builds it.
bool ReadDomain(
    const GU_Detail* gdp,
    bool is2D,
    const int* axes,
    Domain& domain,
    std::string& err);

} // namespace PP

#endif // PP_HOUDINI_DOMAININPUT_H
