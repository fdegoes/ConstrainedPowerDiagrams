#include "houdini/feasibleVolumesSOP.h"
#include "houdini/powerBuilder2DSOP.h"
#include "houdini/powerDiagramSOP.h"
#include "houdini/powerDiagramSolverSOP.h"
#include "houdini/tetrahedralizeSOP.h"
#include "houdini/triangulate2DSOP.h"

#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <UT/UT_DSOVersion.h>

void
newSopOperator(OP_OperatorTable* table)
{
    table->addOperator(new OP_Operator(
        "pp_powerdiagram",
        "PP Power Diagram",
        PP::PowerDiagramSOP::myConstructor,
        PP::PowerDiagramSOP::myTemplateList,
        2,
        2,
        nullptr));

    table->addOperator(new OP_Operator(
        "pp_triangulate2d",
        "PP Triangulate 2D",
        PP::Triangulate2DSOP::myConstructor,
        PP::Triangulate2DSOP::myTemplateList,
        1,
        1,
        nullptr));

    table->addOperator(new OP_Operator(
        "pp_tetrahedralize",
        "PP Tetrahedralize",
        PP::TetrahedralizeSOP::myConstructor,
        PP::TetrahedralizeSOP::myTemplateList,
        1,
        1,
        nullptr));

    table->addOperator(new OP_Operator(
        "pp_powerdiagramsolver",
        "PP Power Diagram Solver",
        PP::PowerDiagramSolverSOP::myConstructor,
        PP::PowerDiagramSolverSOP::GetTemplates(),
        2,
        2,
        nullptr));

    table->addOperator(new OP_Operator(
        "pp_powerbuilder2d",
        "PP Power Builder 2D",
        PP::PowerBuilder2DSOP::myConstructor,
        PP::PowerBuilder2DSOP::myTemplateList,
        2,
        2,
        nullptr));

    table->addOperator(new OP_Operator(
        "pp_feasiblevolumes",
        "PP Feasible Volumes",
        PP::FeasibleVolumesSOP::myConstructor,
        PP::FeasibleVolumesSOP::GetTemplates(),
        2,
        2,
        nullptr));
}
