#ifndef PP_HOUDINI_TRIANGULATE2DSOP_H
#define PP_HOUDINI_TRIANGULATE2DSOP_H

#include <SOP/SOP_Node.h>

namespace PP
{

// Constrained Delaunay triangulation of the region bounded by the
// closed polygons of input 0, lying in the XY, YZ or ZX plane.
//
// Each polygon is a boundary loop; inside/outside follows the even-odd
// rule, so loops inside an outer loop are holes. The output triangles only
// use the polygon points and are the 2D domain of PowerDiagramSOP.
class Triangulate2DSOP : public SOP_Node
{
public:
    static PRM_Template myTemplateList[];
    static OP_Node* myConstructor(OP_Network*, const char*, OP_Operator*);

    Triangulate2DSOP(OP_Network* net, const char* name, OP_Operator* op);
    ~Triangulate2DSOP() override;

protected:
    OP_ERROR cookMySop(OP_Context& context) override;
    const char* inputLabel(OP_InputIdx idx) const override;

private:
    int _PLANE() { return evalInt("plane", 0, 0); }
};

} // namespace PP

#endif // PP_HOUDINI_TRIANGULATE2DSOP_H
