#ifndef PP_HOUDINI_PLANE_H
#define PP_HOUDINI_PLANE_H

#include <PRM/PRM_ChoiceList.h>
#include <PRM/PRM_Name.h>

namespace PP
{

// Axis-aligned plane used by the 2D SOPs ("plane" menu parameter).
enum Plane { PLANE_XY = 0, PLANE_YZ = 1, PLANE_ZX = 2 };

// In-plane axes (a, b) and normal axis n of each plane, so that (a, b, n)
// is right-handed.
inline const int*
GetPlaneAxes(int plane)
{
    static const int axes[3][3] = {
        {0, 1, 2},  // XY
        {1, 2, 0},  // YZ
        {2, 0, 1},  // ZX
    };
    return axes[(plane >= 0 and plane < 3) ? plane : 0];
}

extern PRM_Name       planeParmName;
extern PRM_ChoiceList planeParmMenu;

} // namespace PP

#endif // PP_HOUDINI_PLANE_H
