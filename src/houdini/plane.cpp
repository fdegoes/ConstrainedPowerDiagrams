#include "houdini/plane.h"

namespace PP
{

namespace
{

PRM_Name planeChoices[] = {
    PRM_Name("xy", "XY"),
    PRM_Name("yz", "YZ"),
    PRM_Name("zx", "ZX"),
    PRM_Name(),
};

} // namespace

PRM_Name       planeParmName("plane", "Plane");
PRM_ChoiceList planeParmMenu(PRM_CHOICELIST_SINGLE, planeChoices);

} // namespace PP
