#include "core/lineSearch.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace PP
{

bool
LineSearch::_HasEmptyCell(const ShiftedDiagram& diagram) const
{
    return not (diagram.GetMinVolume() > _params._minVolume);
}

LineSearch::Result
LineSearch::Search(
    const JointDual& dual,
    double value,
    const Vector& g,
    double residual,
    const Vector& dz,
    Vector& z,
    ShiftedDiagram& diagram) const
{
    Result res;
    double omega = 1.0;
    ShiftedDiagram trial;
    Vector gTrial;
    const double roundoff = 16.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::fabs(value));
    const double slope = g.dot(dz);
    for (int k = 0; k <= _params._maxHalvings; ++k, omega *= 0.5) {
        res._halvings = k;
        const Vector zTrial = z + omega * dz;
        if (not dual.Build(zTrial, trial, res._error)) {
            return res;
        }
        if (_HasEmptyCell(trial)) {
            continue;
        }
        const double v = dual.GetValue(zTrial, trial);
        bool accept = v > value;
        if (omega * slope <= roundoff) {
            // The objective change is not measurable: use the residual.
            dual.GetGradient(zTrial, trial, gTrial);
            accept = v >= value - roundoff and dual.GetResidual(gTrial) < residual;
        }
        if (not accept) {
            continue;
        }
        res._accepted = true;
        res._stepLength = omega;
        res._value = v;
        z = zTrial;
        diagram = std::move(trial);
        return res;
    }
    return res;
}

} // namespace PP
