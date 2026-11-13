#ifndef PP_HOUDINI_FEASIBLEVOLUMESSOP_H
#define PP_HOUDINI_FEASIBLEVOLUMESSOP_H

#include "houdini/domainInput.h"

#include <SOP/SOP_Node.h>

#include <vector>

namespace PP
{

// Computes volume fractions for the points of input 0 that are
// consistent with them being the cell centroids of a partition of the
// domain of input 1: f_i > 0, sum f_i = 1, sum f_i P_i = domain centroid
// (the maximum-entropy such fractions, latex/contributions.tex,
// "Remarks / Feasibility").
//
// Works in 2D (triangle mesh in an axis plane) and 3D (tet mesh).
// The fractions are written to a point attribute, ready to be used as the
// targets of PP Power Diagram Solver.
class FeasibleVolumesSOP : public SOP_Node
{
public:
    static PRM_Template* GetTemplates();
    static OP_Node* myConstructor(OP_Network*, const char*, OP_Operator*);

    FeasibleVolumesSOP(OP_Network* net, const char* name, OP_Operator* op);
    ~FeasibleVolumesSOP() override;

protected:
    OP_ERROR cookMySop(OP_Context& context) override;
    bool updateParmsFlags() override;
    const char* inputLabel(OP_InputIdx idx) const override;

private:
    // Writes the fractions to point attribute attrib.
    bool _WriteFractions(const SiteInput& in, const std::vector<double>& fractions, const char* attrib);

    int _MODE() { return evalInt("mode", 0, 0); }
    int _PLANE() { return evalInt("plane", 0, 0); }
    void _OUTATTRIB(UT_String& s, fpreal t) { evalString(s, "outattrib", 0, t); }
    fpreal _TOLERANCE(fpreal t) { return evalFloat("tolerance", 0, t); }
};

} // namespace PP

#endif // PP_HOUDINI_FEASIBLEVOLUMESSOP_H
