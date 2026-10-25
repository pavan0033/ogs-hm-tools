// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>
#include "BaseLib/Error.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
class SourceTermBase
{
public:
    virtual void integrate(const double t, GlobalVector const& x,
                           GlobalVector& b, GlobalMatrix* jac) const = 0;

    /// Commit history only after global acceptance; assembly remains a pure trial.
    virtual bool hasHistory() const { return false; }

    virtual bool supportsSafety() const { return !hasHistory(); }
    virtual void setStrengthReductionFactor(double) {}
    virtual void updateActiveElements(std::vector<std::size_t> const&,bool) {}

    virtual bool supportsCheckpoint() const { return !hasHistory(); }
    virtual std::vector<double> checkpointState() const
    {
        if(hasHistory()) { OGS_FATAL("Structural source does not provide checkpoint state."); }
        return {};
    }
    virtual void restoreCheckpointState(std::vector<double> const& state)
    {
        if(!state.empty() || hasHistory()) { OGS_FATAL("Structural checkpoint state is unsupported."); }
    }

    virtual void postAcceptedTimestep(double, GlobalVector const&) {}

    virtual ~SourceTermBase() = default;
};

class SourceTerm : public SourceTermBase
{
public:
    explicit SourceTerm(
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table)
        : _source_term_dof_table{std::move(source_term_dof_table)}
    {
    }

protected:
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> const _source_term_dof_table;
};

}  // namespace ProcessLib
