// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SourceTermCollection.h"

#include <algorithm>

#include <range/v3/view/filter.hpp>

namespace ProcessLib
{
void SourceTermCollection::addSourceTermsForProcessVariables(
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    unsigned const integration_order,
    const MeshLib::Mesh& bulk_mesh)
{
    for (int variable_id = 0;
         variable_id < static_cast<int>(process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = process_variables[variable_id];
        auto sts =
            pv.createSourceTerms(dof_table, variable_id, integration_order,
                                 _parameters, process_variables, bulk_mesh);

        std::move(sts.begin(), sts.end(), std::back_inserter(_source_terms));
    }
}

void SourceTermCollection::updateActiveElements(std::vector<std::size_t> const& ids,bool masked)
{
    for(auto const& st:_source_terms) { if(st) { st->updateActiveElements(ids,masked); } }
}

bool SourceTermCollection::hasHistory() const
{
    return std::any_of(_source_terms.begin(), _source_terms.end(),
        [](auto const& st) { return st && st->hasHistory(); });
}

bool SourceTermCollection::supportsSafety() const
{
    return std::all_of(_source_terms.begin(),_source_terms.end(),
        [](auto const& st){return !st || st->supportsSafety();});
}
void SourceTermCollection::setSafetyState(bool active,double factor,double physical_time)
{
    safety_active_=active;safety_physical_time_=physical_time;safety_factor_=active?factor:1.;
    for(auto const& st:_source_terms)if(st) { st->setStrengthReductionFactor(safety_factor_); }
}

bool SourceTermCollection::supportsCheckpoint() const
{
    return std::all_of(_source_terms.begin(),_source_terms.end(),
        [](auto const& st){ return !st || st->supportsCheckpoint(); });
}
std::vector<std::vector<double>> SourceTermCollection::checkpointState() const
{
    if(!hasHistory()) { return {}; }
    std::vector<std::vector<double>> result;
    for(auto const& st:_source_terms) { result.push_back(st?st->checkpointState():std::vector<double>{}); }
    return result;
}
void SourceTermCollection::restoreCheckpointState(std::vector<std::vector<double>> const& state)
{
    if(state.empty() && !hasHistory()) { return; }
    if(state.size()!=_source_terms.size()) { OGS_FATAL("Structural checkpoint source count mismatch."); }
    for(std::size_t i=0;i<state.size();++i)
    {
        if(_source_terms[i]) { _source_terms[i]->restoreCheckpointState(state[i]); }
        else if(!state[i].empty()) { OGS_FATAL("Structural checkpoint partition mismatch."); }
    }
}

void SourceTermCollection::postAcceptedTimestep(double const t, GlobalVector const& x)
{
    if(safety_active_ && !_source_terms.empty())
    { INFO("Accepted structural Safety context: coordinate={}, physical_time={}, strength_reduction_factor={}",t,safety_physical_time_,safety_factor_); }
    for (auto const& st : _source_terms)
    {
        if (st) { st->postAcceptedTimestep(safety_active_?safety_physical_time_:t, x); }
    }
}

void SourceTermCollection::integrate(const double t, GlobalVector const& x,
                                     GlobalVector& b, GlobalMatrix* jac) const
{
    // For parallel computing with DDC, a partition may not have source term
    // but a nullptr is assigned to its element in _source_terms.
    auto non_nullptr = [](std::unique_ptr<SourceTermBase> const& st)
    { return st != nullptr; };

    for (auto const& st : _source_terms | ranges::views::filter(non_nullptr))
    {
        st->integrate(safety_active_?safety_physical_time_:t, x, b, jac);
    }
}

}  // namespace ProcessLib
