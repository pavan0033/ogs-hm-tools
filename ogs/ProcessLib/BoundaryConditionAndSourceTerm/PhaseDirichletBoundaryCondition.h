// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <cmath>
#include <memory>
#include "BaseLib/Error.h"
#include "BoundaryCondition.h"
#include "NumLib/IndexValueVector.h"
#include "NumLib/TimeStepping/Time.h"

namespace ProcessLib
{
/// A phase owns (start, end]. At a shared endpoint the preceding phase still
/// owns the accepted state. This avoids imposing the next drainage condition
/// prematurely on the last step of an undrained phase.
class PhaseDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    PhaseDirichletBoundaryCondition(std::unique_ptr<BoundaryCondition> condition,
                                    double const start, double const end)
        : _condition(std::move(condition)), _start(start), _end(end)
    {
        if (!_condition || !std::isfinite(start) || !std::isfinite(end) ||
            NumLib::Time{end} <= NumLib::Time{start})
        {
            OGS_FATAL("PhaseDirichlet requires a condition and finite start_time < end_time.");
        }
    }
    bool requiresExclusiveEssentialDofs() const override { return true; }

    void getEssentialBCValues(
        double const t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& values) const override
    {
        values.ids.clear();
        values.values.clear();
        if (NumLib::Time{t} > NumLib::Time{_start} &&
            NumLib::Time{t} <= NumLib::Time{_end})
        {
            _condition->getEssentialBCValues(t, x, values);
        }
    }
private:
    std::unique_ptr<BoundaryCondition> _condition;
    double const _start;
    double const _end;
};
}  // namespace ProcessLib
