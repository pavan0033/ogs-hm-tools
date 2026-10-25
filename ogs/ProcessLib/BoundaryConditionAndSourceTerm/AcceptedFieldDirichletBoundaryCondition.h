// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include "InitialFieldDirichletBoundaryCondition.h"

namespace ProcessLib
{
/// Hold the accepted primary pressure field during (start,end]. This is a
/// distributed Dirichlet idealisation, not a no-flow or undrained constraint.
class AcceptedFieldDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    AcceptedFieldDirichletBoundaryCondition(
        std::unique_ptr<BoundaryCondition> selector, double const start,
        double const end)
        : _selector(std::move(selector)), _start(start), _end(end)
    {
        if (!_selector || !std::isfinite(start) || !std::isfinite(end) ||
            start<=0 || end<=start)
        { OGS_FATAL("AcceptedFieldDirichlet requires a positive phase interval."); }
    }
    bool requiresExclusiveEssentialDofs() const override { return true; }
    void getEssentialBCValues(double const t, GlobalVector const&,
        NumLib::IndexValueVector<GlobalIndexType>& values) const override
    {
        values.ids.clear();values.values.clear();
        if (NumLib::Time{t}<=NumLib::Time{_start} || NumLib::Time{t}>NumLib::Time{_end}) { return; }
        if (!_captured) { OGS_FATAL("Accepted pressure reference was not captured on the accepted phase boundary."); }
        values=_reference;
    }
    void postAcceptedTimestep(double const t,
        std::vector<GlobalVector*> const& x, int const process_id) override
    {
        if (_captured || NumLib::Time{t}<NumLib::Time{_start}) { return; }
        if (NumLib::Time{t}!=NumLib::Time{_start})
        { OGS_FATAL("Accepted pressure capture time was skipped."); }
        _selector->getEssentialBCValues(t,*x[process_id],_reference);
        for (std::size_t i=0;i<_reference.ids.size();++i)
        {
            _reference.values[i]=x[process_id]->get(_reference.ids[i]);
            if (!std::isfinite(_reference.values[i])) { OGS_FATAL("Nonfinite accepted pressure field."); }
        }
        _captured=true;
    }
    std::vector<double> checkpointState() const override
    {
        std::vector<double> state{_start,_end,_captured?1.:0.};
        for (std::size_t i=0;i<_reference.ids.size();++i)
        { state.push_back(static_cast<double>(_reference.ids[i]));state.push_back(_reference.values[i]); }
        return state;
    }
    void validateCheckpointState(std::vector<double> const& state, double const time) const override
    {
        bool const expected=NumLib::Time{time}>=NumLib::Time{_start};
        if (state.size()<3 || state[0]!=_start || state[1]!=_end ||
            state[2]!=(expected?1.:0.) ||
            state.size()!=(expected?3+2*_indices.size():3))
        { OGS_FATAL("Accepted pressure checkpoint capture/window/layout mismatch."); }
        for (std::size_t i=0;expected && i<_indices.size();++i)
        { if (state[3+2*i]!=static_cast<double>(_indices[i]) || !std::isfinite(state[4+2*i]))
          { OGS_FATAL("Accepted pressure checkpoint DOF/value mismatch."); } }
    }
    void restoreCheckpointState(std::vector<double> const& state, double const time) override
    {
        validateCheckpointState(state,time);_captured=state[2]==1.;
        _reference.ids.clear();_reference.values.clear();
        for (std::size_t i=0;_captured && i<_indices.size();++i)
        { _reference.ids.push_back(_indices[i]);_reference.values.push_back(state[4+2*i]); }
    }
    void setIndices(std::vector<GlobalIndexType> indices) { _indices=std::move(indices); }
private:
    std::unique_ptr<BoundaryCondition> _selector;
    double const _start,_end;
    bool _captured=false;
    std::vector<GlobalIndexType> _indices;
    NumLib::IndexValueVector<GlobalIndexType> _reference;
};
} // namespace ProcessLib
