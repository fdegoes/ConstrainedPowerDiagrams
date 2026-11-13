#ifndef PP_HOUDINI_POWERDIAGRAMSOP_H
#define PP_HOUDINI_POWERDIAGRAMSOP_H

#include <SOP/SOP_Node.h>

namespace PP
{

// Power diagram of the points of input 0 clipped to the polytope of
// input 1.
//
// Input 0: points P, optional float weight (w) and vector shift (u).
//   The diagram uses sites P + u with weights w + |u|^2.
// Input 1: 2D mode: a triangle mesh lying in the XY, YZ or ZX plane.
//   3D mode: a tetrahedral mesh (tetrahedron primitives).
//
// Output: polygons with a primitive int attribute "cell" holding the index
// of the input point that generated them.
class PowerDiagramSOP : public SOP_Node
{
public:
    static PRM_Template myTemplateList[];
    static OP_Node* myConstructor(OP_Network*, const char*, OP_Operator*);

    PowerDiagramSOP(OP_Network* net, const char* name, OP_Operator* op);
    ~PowerDiagramSOP() override;

protected:
    OP_ERROR cookMySop(OP_Context& context) override;
    bool updateParmsFlags() override;
    const char* inputLabel(OP_InputIdx idx) const override;

private:
    int _MODE() { return evalInt("mode", 0, 0); }
    int _PLANE() { return evalInt("plane", 0, 0); }
    void _WEIGHTATTRIB(UT_String& s, fpreal t) { evalString(s, "weightattrib", 0, t); }
    void _SHIFTATTRIB(UT_String& s, fpreal t) { evalString(s, "shiftattrib", 0, t); }
    fpreal _BOUNDARYANGLE(fpreal t) { return evalFloat("boundaryangle", 0, t); }
};

} // namespace PP

#endif // PP_HOUDINI_POWERDIAGRAMSOP_H
