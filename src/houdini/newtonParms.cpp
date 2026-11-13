#include "houdini/newtonParms.h"

#include <OP/OP_Parameters.h>
#include <PRM/PRM_Include.h>

#include <iostream>

namespace PP
{

namespace
{

PRM_Name maxIterName("maxiterations", "Max Iterations");
PRM_Default maxIterDefault(100);
PRM_Range maxIterRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 500);

PRM_Name regName("regularization", "Regularization");
PRM_Default regDefault(1e-3);
PRM_Range regRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_UI, 1e-2);

PRM_Name maxHalvingsName("maxhalvings", "Max Step Halvings");
PRM_Default maxHalvingsDefault(40);
PRM_Range maxHalvingsRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 100);

PRM_Name logName("verbose", "Print Iterations");

PRM_Name separatorName("newtonsep", "");

} // namespace

void
AppendNewtonParmTemplates(std::vector<PRM_Template>& templates)
{
    templates.push_back(PRM_Template(PRM_INT, 1, &maxIterName, &maxIterDefault, nullptr, &maxIterRange));
    templates.push_back(PRM_Template(PRM_FLT_LOG, 1, &regName, &regDefault, nullptr, &regRange));
    templates.push_back(PRM_Template(PRM_INT, 1, &maxHalvingsName, &maxHalvingsDefault, nullptr,
                                     &maxHalvingsRange));
    templates.push_back(PRM_Template(PRM_SEPARATOR, 1, &separatorName));
    templates.push_back(PRM_Template(PRM_TOGGLE, 1, &logName, PRMzeroDefaults));
}

NewtonParms
ReadNewtonParameters(OP_Parameters* node, fpreal t)
{
    NewtonParms p;
    p._maxIterations = int(node->evalInt("maxiterations", 0, t));
    p._regularization = node->evalFloat("regularization", 0, t);
    p._maxHalvings = int(node->evalInt("maxhalvings", 0, t));
    p._log = node->evalInt("verbose", 0, t) ? &std::cout : nullptr;
    return p;
}

} // namespace PP
