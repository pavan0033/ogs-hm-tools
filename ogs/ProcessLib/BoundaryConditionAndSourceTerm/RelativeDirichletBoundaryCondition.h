// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <algorithm>
#include <cmath>
#include <memory>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "BoundaryCondition.h"
#include "NumLib/IndexValueVector.h"

namespace ProcessLib
{
/// Inactive through reference_time; afterwards x = accepted_reference + increment.
/// Neither essential-value queries nor rejected trials may change the reference.
class RelativeDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    RelativeDirichletBoundaryCondition(
        std::unique_ptr<BoundaryCondition> increment,
        double const reference_time)
        : _increment(std::move(increment)), _reference_time(reference_time)
    {
        if (!_increment || !std::isfinite(reference_time) || reference_time <= 0)
        {
            OGS_FATAL("RelativeDirichlet requires finite positive reference_time.");
        }
    }

    bool requiresExclusiveEssentialDofs() const override { return true; }

    void getEssentialBCValues(
        double const t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& values) const override
    {
        values.ids.clear();
        values.values.clear();
        if (t <= _reference_time + timeTolerance())
        {
            return;
        }
        if (!_captured)
        {
            OGS_FATAL("RelativeDirichlet reference was not captured at an accepted "
                      "reference_time. Include that time in the time grid.");
        }
        _increment->getEssentialBCValues(t, x, values);
        if (values.ids != _reference.ids)
        {
            OGS_FATAL("RelativeDirichlet DOF map changed after capture.");
        }
        for (std::size_t i = 0; i < values.values.size(); ++i)
        {
            values.values[i] += _reference.values[i];
            if (!std::isfinite(values.values[i]))
            {
                OGS_FATAL("RelativeDirichlet target is nonfinite.");
            }
        }
    }

    void postAcceptedTimestep(double const t,
                              std::vector<GlobalVector*> const& x,
                              int const process_id) override
    {
        if (_captured || t < _reference_time - timeTolerance())
        {
            return;
        }
        if (std::abs(t - _reference_time) > timeTolerance())
        {
            OGS_FATAL("RelativeDirichlet reference_time was skipped by the "
                      "accepted time grid.");
        }
        auto const& solution = *x[process_id];
        _increment->getEssentialBCValues(t, solution, _reference);
        for (std::size_t i = 0; i < _reference.ids.size(); ++i)
        {
            if (!std::isfinite(_reference.values[i]) ||
                std::abs(_reference.values[i]) > 1e-14)
            {
                OGS_FATAL("RelativeDirichlet increment must be zero at "
                          "reference_time, in displacement units.");
            }
            _reference.values[i] = solution.get(_reference.ids[i]);
            if (!std::isfinite(_reference.values[i]))
            {
                OGS_FATAL("RelativeDirichlet accepted reference is nonfinite.");
            }
        }
        _captured = true;
        INFO("RelativeDirichlet captured accepted reference at t={} for {} DOFs.",
             t, _reference.ids.size());
    }

private:
    double timeTolerance() const
    {
        return 1e-12 * std::max(1., std::abs(_reference_time));
    }
    std::unique_ptr<BoundaryCondition> _increment;
    double const _reference_time;
    bool _captured = false;
    NumLib::IndexValueVector<GlobalIndexType> _reference;
};
}  // namespace ProcessLib
