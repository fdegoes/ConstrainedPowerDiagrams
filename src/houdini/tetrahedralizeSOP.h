#ifndef PP_HOUDINI_TETRAHEDRALIZESOP_H
#define PP_HOUDINI_TETRAHEDRALIZESOP_H

#include "core/tetrahedralize.h"

#include <SOP/SOP_Node.h>

namespace PP
{

// Fills the volume bounded by the triangulated surface of input 0
// with tetrahedra (TetGen through geogram).
//
// The output tetrahedron primitives are the 3D domain of PowerDiagramSOP.
class TetrahedralizeSOP : public SOP_Node
{
public:
    static PRM_Template myTemplateList[];
    static OP_Node* myConstructor(OP_Network*, const char*, OP_Operator*);

    TetrahedralizeSOP(OP_Network* net, const char* name, OP_Operator* op);
    ~TetrahedralizeSOP() override;

protected:
    OP_ERROR cookMySop(OP_Context& context) override;
    bool updateParmsFlags() override;
    const char* inputLabel(OP_InputIdx idx) const override;

private:
    TetrahedralizeOptions _ReadOptions(fpreal t);

    int _PREPROCESS() { return evalInt("preprocess", 0, 0); }
    int _REFINE() { return evalInt("refine", 0, 0); }
    fpreal _QUALITY(fpreal t) { return evalFloat("quality", 0, t); }
};

} // namespace PP

#endif // PP_HOUDINI_TETRAHEDRALIZESOP_H
