#ifndef PP_HOUDINI_POWERDIAGRAMSOLVERSOP_H
#define PP_HOUDINI_POWERDIAGRAMSOLVERSOP_H

#include "core/dualSolver.h"
#include "houdini/domainInput.h"

#include <SOP/SOP_Node.h>

#include <vector>

namespace PP
{

// Solves for the weights of a power diagram whose cells, clipped to
// the domain of input 1, satisfy prescribed constraints.
//
// Constraint "Volume": each site gets a target volume fraction (point
// attribute, fractions summing to 1); the target cell volumes are the
// fractions times the domain volume. Newton's method maximizes the concave
// dual energy of latex/background.tex ("Volume constraints").
//
// Output: input 0 with the solved weights. Sites are P + u and the solver
// works on the effective weight w + |u|^2, so feeding the output to
// PP Power Diagram (same shift attribute) reproduces the solved cells.
class PowerDiagramSolverSOP : public SOP_Node
{
public:
    static PRM_Template* GetTemplates();
    static OP_Node* myConstructor(OP_Network*, const char*, OP_Operator*);

    PowerDiagramSolverSOP(OP_Network* net, const char* name, OP_Operator* op);
    ~PowerDiagramSolverSOP() override;

protected:
    OP_ERROR cookMySop(OP_Context& context) override;
    bool updateParmsFlags() override;
    const char* inputLabel(OP_InputIdx idx) const override;

private:
    enum Constraint { CONSTRAINT_VOLUME = 0 };

    // Target fractions of attribute attrib, normalized and scaled to
    // volumes. Adds an error and returns false if they are missing or not
    // positive.
    bool _ReadTargetVolumes(
        const SiteInput& in,
        const char* attrib,
        double domainVolume,
        std::vector<double>& targets);

    DualSolver::Parameters _ReadSolverParameters(fpreal t);

    // Warnings for a non-converged solve.
    void _ReportSolve(const DualSolver::Result& res);

    // Point attributes (weights, achieved fractions) and detail attributes
    // (solver status).
    bool _WriteWeights(
        const SiteInput& in,
        const std::vector<double>& omega,
        const DualSolver::Result& res,
        const char* attrib,
        double domainVolume);

    int _MODE() { return evalInt("mode", 0, 0); }
    int _PLANE() { return evalInt("plane", 0, 0); }
    int _CONSTRAINT() { return evalInt("constraint", 0, 0); }
    void _TARGETATTRIB(UT_String& s, fpreal t) { evalString(s, "targetattrib", 0, t); }
    void _WEIGHTATTRIB(UT_String& s, fpreal t) { evalString(s, "weightattrib", 0, t); }
    void _SHIFTATTRIB(UT_String& s, fpreal t) { evalString(s, "shiftattrib", 0, t); }
    void _OUTWEIGHTATTRIB(UT_String& s, fpreal t) { evalString(s, "outweightattrib", 0, t); }
    fpreal _VOLUMETOL(fpreal t) { return evalFloat("volumetolerance", 0, t); }
};

} // namespace PP

#endif // PP_HOUDINI_POWERDIAGRAMSOLVERSOP_H
