#ifndef PP_HOUDINI_POWERBUILDER2DSOP_H
#define PP_HOUDINI_POWERBUILDER2DSOP_H

#include <SOP/SOP_Node.h>

namespace PP
{

// Power diagram of shifted weighted points, clipped to a 2D triangle mesh.
//
// Input 0: points with P = q, and optional point attributes for the weights
//   w (float) and shifts u (vector), named by parameters; an empty or
//   missing attribute counts as 0. The generators are p = q + u (in-plane
//   part of u) with weights w: the power distance of x to point i is
//   |x - p_i|^2 - w_i.
// Input 1: a triangle mesh in the XY, YZ or ZX plane, the domain.
//
// Output: one closed polygon per non-empty cell, on its own points (points
// shared between cells are duplicated), with a primitive int attribute
// "cell" (index of the input point). Optional primitive attributes, named by
// parameters: the cell volume (area, float) and centroid (vector, in the
// plane at the height of the domain). Cells are assumed connected and
// without holes (one boundary loop each); otherwise the cook fails.
class PowerBuilder2DSOP : public SOP_Node
{
public:
    static PRM_Template myTemplateList[];
    static OP_Node* myConstructor(OP_Network*, const char*, OP_Operator*);

    PowerBuilder2DSOP(OP_Network* net, const char* name, OP_Operator* op);
    ~PowerBuilder2DSOP() override;

protected:
    OP_ERROR cookMySop(OP_Context& context) override;
    const char* inputLabel(OP_InputIdx idx) const override;

private:
    // Generators p = q + u and weights w of the points of gdp (2D, along
    // axes). Adds warnings for missing attributes.
    void _ReadGenerators(
        const GU_Detail* gdp,
        const int* axes,
        const char* weightAttrib,
        const char* shiftAttrib,
        std::vector<double>& generators,
        std::vector<double>& weights);

    // Primitive attributes of the cells: "cell", and the optional volume and
    // centroid (empty names skip them).
    void _WriteCellAttributes(
        GA_Offset firstPrim,
        const std::vector<double>& pts,
        const std::vector<int>& polySizes,
        const std::vector<int>& polyVerts,
        const std::vector<int>& polyCell,
        const int* axes,
        double height,
        const char* volumeAttrib,
        const char* centroidAttrib);

    int _PLANE() { return evalInt("plane", 0, 0); }
    void _WEIGHTATTRIB(UT_String& s, fpreal t) { evalString(s, "weightattrib", 0, t); }
    void _SHIFTATTRIB(UT_String& s, fpreal t) { evalString(s, "shiftattrib", 0, t); }
    void _VOLUMEATTRIB(UT_String& s, fpreal t) { evalString(s, "volumeattrib", 0, t); }
    void _CENTROIDATTRIB(UT_String& s, fpreal t) { evalString(s, "centroidattrib", 0, t); }
};

} // namespace PP

#endif // PP_HOUDINI_POWERBUILDER2DSOP_H
