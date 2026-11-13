#include "houdini/powerBuilder2DSOP.h"
#include "houdini/domainInput.h"
#include "houdini/geometryIO.h"
#include "houdini/plane.h"

#include "core/powerDiagram.h"

#include <GA/GA_Handle.h>
#include <GU/GU_Detail.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <PRM/PRM_Include.h>

#include <string>
#include <vector>

namespace PP
{

namespace
{

PRM_Name weightAttribName("weightattrib", "Weight Attribute");
PRM_Default weightAttribDefault(0, "w");

PRM_Name shiftAttribName("shiftattrib", "Shift Attribute");
PRM_Default shiftAttribDefault(0, "u");

PRM_Name volumeAttribName("volumeattrib", "Volume Attribute");
PRM_Default volumeAttribDefault(0, "volume");

PRM_Name centroidAttribName("centroidattrib", "Centroid Attribute");
PRM_Default centroidAttribDefault(0, "centroid");

PRM_Name separatorName("outputsep", "");

// Area and centroid of the closed 2D polygon with point indices
// verts[first, first + size) into pts (shoelace formula).
void
_ComputePolygonMoments(
    const std::vector<double>& pts,
    const std::vector<int>& verts,
    size_t first,
    int size,
    double& area,
    double centroid[2])
{
    double a = 0.0, cx = 0.0, cy = 0.0;
    for (int k = 0; k < size; ++k) {
        const double* p = &pts[2 * verts[first + k]];
        const double* q = &pts[2 * verts[first + (k + 1) % size]];
        const double cross = p[0] * q[1] - q[0] * p[1];
        a += cross;
        cx += (p[0] + q[0]) * cross;
        cy += (p[1] + q[1]) * cross;
    }
    area = 0.5 * a;
    centroid[0] = area != 0.0 ? cx / (6.0 * area) : 0.0;
    centroid[1] = area != 0.0 ? cy / (6.0 * area) : 0.0;
}

// True if some cell owns several polygons (several pieces or holes).
bool
_HasMultiLoopCell(const std::vector<int>& polyCell, size_t nbCells)
{
    std::vector<char> seen(nbCells, 0);
    for (int c : polyCell) {
        if (seen[c]) {
            return true;
        }
        seen[c] = 1;
    }
    return false;
}

} // namespace

PRM_Template PowerBuilder2DSOP::myTemplateList[] = {
    PRM_Template(PRM_ORD, 1, &planeParmName, PRMzeroDefaults, &planeParmMenu),
    PRM_Template(PRM_STRING, 1, &weightAttribName, &weightAttribDefault),
    PRM_Template(PRM_STRING, 1, &shiftAttribName, &shiftAttribDefault),
    PRM_Template(PRM_SEPARATOR, 1, &separatorName),
    PRM_Template(PRM_STRING, 1, &volumeAttribName, &volumeAttribDefault),
    PRM_Template(PRM_STRING, 1, &centroidAttribName, &centroidAttribDefault),
    PRM_Template(),
};

OP_Node*
PowerBuilder2DSOP::myConstructor(OP_Network* net, const char* name, OP_Operator* op)
{
    return new PowerBuilder2DSOP(net, name, op);
}

PowerBuilder2DSOP::PowerBuilder2DSOP(OP_Network* net, const char* name, OP_Operator* op)
    : SOP_Node(net, name, op)
{ }

PowerBuilder2DSOP::~PowerBuilder2DSOP() { }

const char*
PowerBuilder2DSOP::inputLabel(OP_InputIdx idx) const
{
    switch (idx) {
    case 0: return "Points (P = q, w, u)";
    case 1: return "Domain (triangles)";
    default: return "Invalid";
    }
}

void
PowerBuilder2DSOP::_ReadGenerators(
    const GU_Detail* gdp,
    const int* axes,
    const char* weightAttrib,
    const char* shiftAttrib,
    std::vector<double>& generators,
    std::vector<double>& weights)
{
    GA_ROHandleF weightH;
    if (*weightAttrib) {
        weightH = GA_ROHandleF(gdp, GA_ATTRIB_POINT, weightAttrib);
        if (not weightH.isValid()) {
            addWarning(SOP_MESSAGE, (std::string("No float point attribute '") + weightAttrib +
                                     "'; the weights are 0.").c_str());
        }
    }
    GA_ROHandleV3 shiftH;
    if (*shiftAttrib) {
        shiftH = GA_ROHandleV3(gdp, GA_ATTRIB_POINT, shiftAttrib);
        if (not shiftH.isValid()) {
            addWarning(SOP_MESSAGE, (std::string("No vector point attribute '") + shiftAttrib +
                                     "'; the shifts are 0.").c_str());
        }
    }
    generators.clear();
    weights.clear();
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff)
    {
        UT_Vector3D p(gdp->getPos3(ptoff));
        if (shiftH.isValid()) {
            p += UT_Vector3D(shiftH.get(ptoff));
        }
        generators.push_back(p[axes[0]]);
        generators.push_back(p[axes[1]]);
        weights.push_back(weightH.isValid() ? double(weightH.get(ptoff)) : 0.0);
    }
}

void
PowerBuilder2DSOP::_WriteCellAttributes(
    GA_Offset firstPrim,
    const std::vector<double>& pts,
    const std::vector<int>& polySizes,
    const std::vector<int>& polyVerts,
    const std::vector<int>& polyCell,
    const int* axes,
    double height,
    const char* volumeAttrib,
    const char* centroidAttrib)
{
    GA_RWHandleI cellH(gdp->addIntTuple(GA_ATTRIB_PRIMITIVE, "cell", 1));
    GA_RWHandleF volumeH;
    if (*volumeAttrib) {
        volumeH = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, volumeAttrib, 1));
    }
    GA_RWHandleV3 centroidH;
    if (*centroidAttrib) {
        centroidH = GA_RWHandleV3(gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, centroidAttrib, 3));
        if (centroidH.isValid()) {
            centroidH->setTypeInfo(GA_TYPE_POINT);
        }
    }
    size_t first = 0;
    for (size_t k = 0; k < polySizes.size(); ++k) {
        const GA_Offset prim = firstPrim + GA_Offset(k);
        cellH.set(prim, polyCell[k]);
        double area = 0.0, c[2] = {0.0, 0.0};
        _ComputePolygonMoments(pts, polyVerts, first, polySizes[k], area, c);
        if (volumeH.isValid()) {
            volumeH.set(prim, float(area));
        }
        if (centroidH.isValid()) {
            UT_Vector3 pos;
            pos[axes[0]] = float(c[0]);
            pos[axes[1]] = float(c[1]);
            pos[axes[2]] = float(height);
            centroidH.set(prim, pos);
        }
        first += size_t(polySizes[k]);
    }
}

OP_ERROR
PowerBuilder2DSOP::cookMySop(OP_Context& context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT) {
        return error();
    }

    const fpreal t = context.getTime();
    gdp->clearAndDestroy();
    const int* axes = GetPlaneAxes(_PLANE());

    UT_String weightName, shiftName, volumeName, centroidName;
    _WEIGHTATTRIB(weightName, t);
    _SHIFTATTRIB(shiftName, t);
    _VOLUMEATTRIB(volumeName, t);
    _CENTROIDATTRIB(centroidName, t);

    std::vector<double> generators, weights;
    _ReadGenerators(inputGeo(0, context), axes, weightName.c_str(), shiftName.c_str(), generators, weights);
    if (weights.empty()) {
        addError(SOP_MESSAGE, "The first input has no points.");
        return error();
    }

    std::vector<double> domPts;
    std::vector<int> tris;
    double height = 0.0;
    std::string err;
    if (not ReadDomain(inputGeo(1, context), true, axes, domPts, tris, height, err)) {
        addError(SOP_MESSAGE, err.c_str());
        return error();
    }

    PolyMesh cells;
    if (not ComputePowerDiagram2D(generators, weights, domPts, tris, cells, err)) {
        addError(SOP_MESSAGE, err.c_str());
        return error();
    }
    if (_HasMultiLoopCell(cells._polyCell, weights.size())) {
        addError(SOP_MESSAGE, "A cell has several pieces or a hole; Power Builder 2D assumes one "
                              "boundary loop per cell (use PP Power Diagram for such domains).");
        return error();
    }

    const GA_Offset firstPrim =
        AppendPolygons(gdp, MakePositions3D(cells._points, 2, axes, height), cells._polySizes, cells._polyVerts);
    _WriteCellAttributes(firstPrim, cells._points, cells._polySizes, cells._polyVerts, cells._polyCell, axes,
                         height, volumeName.c_str(), centroidName.c_str());
    SetDetailInt(gdp, "nb_empty_cells", int(weights.size() - cells._polySizes.size()));
    return error();
}

} // namespace PP
