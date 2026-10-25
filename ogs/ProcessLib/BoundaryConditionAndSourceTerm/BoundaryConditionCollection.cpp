// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BoundaryConditionCollection.h"

#include <typeinfo>

#include "MathLib/LinAlg/LinAlg.h"
#include "ReleaseNodalForce.h"
#include "PumpingWellBoundaryCondition.h"
#include "ReservoirSeepageBoundaryCondition.h"

namespace ProcessLib
{

bool BoundaryConditionCollection::hasEnabledActiveSet(double const t) const
{
    for (auto const& bc : _boundary_conditions)
    {
        if (!bc->hasActiveSet())
        {
            continue;
        }
        if (auto const* well =
                dynamic_cast<PumpingWellBoundaryCondition const*>(bc.get()))
        {
            if (well->enabled(t))
            {
                return true;
            }
            continue;
        }
        if (auto const* face =
                dynamic_cast<ReservoirSeepageBoundaryCondition const*>(bc.get()))
        {
            if (face->enabled(t))
            {
                return true;
            }
            continue;
        }
        return true;
    }
    return false;
}

bool BoundaryConditionCollection::updateActiveSet(double t, GlobalVector const& x) const
{
    // HM MassFlowRate is the bulk residual, before natural boundary terms.
    // A well sink must not be counted as atmospheric surface discharge.
    std::map<std::size_t,double> pumping;
    for (auto const& bc : _boundary_conditions)
    {
        if (auto const* well=dynamic_cast<PumpingWellBoundaryCondition const*>(bc.get()))
        { for(auto const& [node,rate] : well->massRhs(t)) { pumping[node]+=rate; } }
    }
    for (auto const& bc : _boundary_conditions)
    {
        if (auto* face=dynamic_cast<ReservoirSeepageBoundaryCondition*>(bc.get()))
        { face->setPumpingMassRhs(pumping); }
    }
    bool changed=false;
    for (auto const& bc : _boundary_conditions)
    { changed=bc->updateActiveSet(t,x)||changed; }
    return changed;
}

void BoundaryConditionCollection::setReleaseNodalForces(
    GlobalVector const* r_neq) const
{
    MathLib::LinAlg::setLocalAccessibleVector(*r_neq);

    for (auto const& bc : _boundary_conditions)
    {
        auto* release_nodal_forces = dynamic_cast<ReleaseNodalForce*>(bc.get());
        if (!release_nodal_forces)
        {
            continue;
        }
        // For ReleasedNodalForce, we need to set the non-equilibrium
        // initial residuum vector.
        release_nodal_forces->set(r_neq);
    }
}

void BoundaryConditionCollection::applyNaturalBC(
    const double t, std::vector<GlobalVector*> const& x, int const process_id,
    GlobalMatrix* K, GlobalVector& b, GlobalMatrix* Jac) const
{
    for (auto const& bc : _boundary_conditions)
    {
        bc->applyNaturalBC(t, x, process_id, K, b, Jac);
    }
}

void BoundaryConditionCollection::addBCsForProcessVariables(
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    unsigned const integration_order, Process const& process,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    for (int variable_id = 0;
         variable_id < static_cast<int>(process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = process_variables[variable_id];
        auto bcs = pv.createBoundaryConditions(
            dof_table, variable_id, integration_order, _parameters, process,
            process_variables, media);

        std::move(bcs.begin(), bcs.end(),
                  std::back_inserter(_boundary_conditions));
    }

    for (std::size_t i=0;i<_boundary_conditions.size();++i)
    {
        auto const* a=dynamic_cast<PumpingWellBoundaryCondition const*>(_boundary_conditions[i].get());
        if (!a) { continue; }
        for (std::size_t j=i+1;j<_boundary_conditions.size();++j)
        {
            auto const* b=dynamic_cast<PumpingWellBoundaryCondition const*>(_boundary_conditions[j].get());
            if (b && a->overlaps(*b))
            { OGS_FATAL("PumpingWell screens cannot share pressure nodes in overlapping time intervals."); }
        }
    }
    // For each BC there will be storage for Dirichlet BC. This storage will be
    // uninitialized by default, and has to be filled by the respective BC
    // object if needed.
    _dirichlet_bcs.resize(_boundary_conditions.size());
}
}  // namespace ProcessLib
