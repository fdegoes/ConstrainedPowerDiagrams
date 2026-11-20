#ifndef PP_CORE_LINESEARCH_H
#define PP_CORE_LINESEARCH_H

#include "core/jointDual.h"
#include "core/shiftedDiagram.h"

#include <string>

namespace PP
{

// Backtracking line search on a concave dual (latex/appendix-algorithms.tex,
// Algorithm "line-search"): from z along dz, halves the step length omega
// (starting at 1) until no cell of the power diagram at z + omega dz is
// empty and the objective increases.
//
// Close to the maximizer the increase of the objective drops below its
// round-off and can no longer be measured. When the predicted increase
// omega g^T dz is within round-off of the objective, the trial is accepted
// instead if it decreases the gradient residual (JointDual::GetResidual).
class LineSearch
{
public:
    using Vector = JointDual::Vector;

    struct Parameters
    {
        int    _maxHalvings = 40;
        // A cell counts as empty when its volume is at most minVolume
        // (0 keeps any non-empty cell).
        double _minVolume = 0.0;
    };

    struct Result
    {
        bool        _accepted = false;
        double      _stepLength = 0.0;  // omega
        int         _halvings = 0;
        double      _value = 0.0;       // objective at the accepted point
        std::string _error;             // why the last unbuilt trial diagram failed, if any
    };

    LineSearch() = default;
    explicit LineSearch(const Parameters& params) : _params(params) { }

    Parameters& GetParameters() { return _params; }
    const Parameters& GetParameters() const { return _params; }

    // value, g and residual are the objective, gradient and residual at z.
    // On success, z and diagram hold the accepted point and its power
    // diagram; otherwise they are unchanged.
    Result Search(
        const JointDual& dual,
        double value,
        const Vector& g,
        double residual,
        const Vector& dz,
        Vector& z,
        ShiftedDiagram& diagram) const;

private:
    bool _HasEmptyCell(const ShiftedDiagram& diagram) const;

    Parameters _params;
};

} // namespace PP

#endif // PP_CORE_LINESEARCH_H
