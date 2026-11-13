#ifndef PP_HOUDINI_GEOMETRYIO_H
#define PP_HOUDINI_GEOMETRYIO_H

// Conversions between Houdini geometry and the raw arrays of the core
// library.

#include <GA/GA_Types.h>
#include <SYS/SYS_Types.h>
#include <UT/UT_Vector3.h>

#include <string>
#include <vector>

class GU_Detail;

namespace PP
{

// Positions of all points of gdp, dim coordinates each, taken along
// axes (the first dim entries).
std::vector<double> ReadPositions(const GU_Detail* gdp, const int* axes, int dim);

// Point indices of all triangles of gdp. Fails if a polygon is not a
// triangle; hint is appended to the error.
bool ReadTriangles(const GU_Detail* gdp, std::vector<int>& tris, const char* hint, std::string& err);

// Point indices of all tetrahedra of gdp (other primitives ignored).
std::vector<int> ReadTetrahedra(const GU_Detail* gdp);

// Average coordinate of the points of all polygons along axis.
double ComputeAverageHeight(const GU_Detail* gdp, int axis);

// Polygons of gdp with at least 3 vertices, as loops of 2D points along
// axes. height is the average coordinate along axes[2].
void ReadLoops(
    const GU_Detail* gdp,
    const int* axes,
    std::vector<double>& loopPts,
    std::vector<int>& loopSizes,
    double& height);

// Converts dim-dimensional points to 3D: in 2D the coordinates go along
// axes[0], axes[1] and axes[2] is set to height.
std::vector<UT_Vector3> MakePositions3D(const std::vector<double>& pts, int dim, const int* axes, double height);

// Appends closed polygons; returns the offset of the first one.
GA_Offset AppendPolygons(
    GU_Detail* gdp,
    const std::vector<UT_Vector3>& pos,
    const std::vector<int>& polySizes,
    const std::vector<int>& polyVerts);

// Appends tetrahedra on new points.
void AppendTetrahedra(GU_Detail* gdp, const std::vector<double>& pts, const std::vector<int>& tets);

void SetDetailInt(GU_Detail* gdp, const char* name, int value);
void SetDetailFloat(GU_Detail* gdp, const char* name, double value);
void SetDetailString(GU_Detail* gdp, const char* name, const char* value);

} // namespace PP

#endif // PP_HOUDINI_GEOMETRYIO_H
