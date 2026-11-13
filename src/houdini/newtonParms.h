#ifndef PP_HOUDINI_NEWTONPARMS_H
#define PP_HOUDINI_NEWTONPARMS_H

#include <PRM/PRM_Template.h>
#include <SYS/SYS_Types.h>

#include <iosfwd>
#include <vector>

class OP_Parameters;

namespace PP
{

// Number of parameters in the "Newton Solver" tab (for PRM_SWITCHER).
constexpr int kNumNewtonParms = 5;

// Values of the "Newton Solver" tab, shared by the solver SOPs.
struct NewtonParms
{
    int           _maxIterations = 100;
    double        _regularization = 1e-3;
    int           _maxHalvings = 40;
    std::ostream* _log = nullptr;  // "Print Iterations" logs to std::cout
};

// Appends the parameters of the "Newton Solver" tab, in order, to
// templates (right after the parameters of the previous tab).
// Call it at operator registration time, not during static initialization.
void AppendNewtonParmTemplates(std::vector<PRM_Template>& templates);

// Reads the "Newton Solver" tab of node.
NewtonParms ReadNewtonParameters(OP_Parameters* node, fpreal t);

} // namespace PP

#endif // PP_HOUDINI_NEWTONPARMS_H
