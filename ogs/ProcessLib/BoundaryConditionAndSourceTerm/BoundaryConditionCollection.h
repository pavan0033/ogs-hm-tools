// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <unordered_set>

#include "BaseLib/Error.h"
#include "BoundaryCondition.h"
#include "MaterialLib/MPL/Medium.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/ProcessVariable.h"

namespace ProcessLib
{
class BoundaryConditionCollection final
{
public:
    explicit BoundaryConditionCollection(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters)
        : _parameters(parameters)
    {
    }

    bool hasActiveSet() const
    {
        for (auto const& bc:_boundary_conditions) { if (bc->hasActiveSet()) { return true; } }
        return false;
    }
    /// Whether a declared active-set boundary needs a fresh residual at t.
    /// Unknown boundary types conservatively retain their existing behaviour.
    bool hasEnabledActiveSet(double t) const;
    bool updateActiveSet(double t, GlobalVector const& x) const;
    void setReleaseNodalForces(GlobalVector const* r_neq) const;

    void applyNaturalBC(const double t, std::vector<GlobalVector*> const& x,
                        int const process_id, GlobalMatrix* K, GlobalVector& b,
                        GlobalMatrix* Jac) const;

    std::vector<NumLib::IndexValueVector<GlobalIndexType>> const*
    getKnownSolutions(double const t, GlobalVector const& x) const
    {
        auto const n_bcs = _boundary_conditions.size();
        for (std::size_t i = 0; i < n_bcs; ++i)
        {
            auto const& bc = *_boundary_conditions[i];
            auto& dirichlet_storage = _dirichlet_bcs[i];
            bc.getEssentialBCValues(t, x, dirichlet_storage);
        }
        for (std::size_t i = 0; i < n_bcs; ++i)
        {
            if (!_boundary_conditions[i]->requiresExclusiveEssentialDofs())
            {
                continue;
            }
            auto const& ids = _dirichlet_bcs[i].ids;
            std::unordered_set<GlobalIndexType> const controlled(ids.begin(), ids.end());
            for (std::size_t j = 0; j < n_bcs; ++j)
            {
                if (j == i)
                {
                    continue;
                }
                for (auto const id : _dirichlet_bcs[j].ids)
                {
                    if (controlled.contains(id))
                    {
                        OGS_FATAL("Phase-controlled boundary has competing "
                                  "essential boundary conditions on DOF {}.", id);
                    }
                }
            }
        }
        return &_dirichlet_bcs;
    }

    void addBCsForProcessVariables(
        std::vector<std::reference_wrapper<ProcessVariable>> const&
            process_variables,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        unsigned const integration_order, Process const& process,
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    void addBoundaryCondition(std::unique_ptr<BoundaryCondition>&& bc)
    {
        _boundary_conditions.push_back(std::move(bc));
    }

    void preTimestep(const double t, std::vector<GlobalVector*> const& x,
                     int const process_id) const
    {
        for (auto const& bc_ptr : _boundary_conditions)
        {
            bc_ptr->preTimestep(t, x, process_id);
        }
    }

    void postTimestep(const double t, std::vector<GlobalVector*> const& x,
                      int const process_id) const
    {
        for (auto const& bc_ptr : _boundary_conditions)
        {
            bc_ptr->postTimestep(t, x, process_id);
        }
    }

    void postAcceptedTimestep(double const t,
                              std::vector<GlobalVector*> const& x,
                              int const process_id) const
    {
        for (auto const& bc : _boundary_conditions)
        {
            bc->postAcceptedTimestep(t, x, process_id);
        }
    }

    std::vector<std::vector<double>> checkpointState() const
    {
        std::vector<std::vector<double>> states;
        for (auto const& bc : _boundary_conditions) { states.push_back(bc->checkpointState()); }
        return states;
    }
    void validateCheckpointState(std::vector<std::vector<double>> const& states,
                                  double const time) const
    {
        if (states.size() != _boundary_conditions.size())
        { OGS_FATAL("Boundary checkpoint catalogue size mismatch."); }
        for (std::size_t i=0; i<states.size(); ++i)
        { _boundary_conditions[i]->validateCheckpointState(states[i],time); }
    }
    void restoreCheckpointState(std::vector<std::vector<double>> const& states,
                                 double const time) const
    {
        validateCheckpointState(states,time);
        for (std::size_t i=0; i<states.size(); ++i)
        { _boundary_conditions[i]->restoreCheckpointState(states[i],time); }
    }

private:
    mutable std::vector<NumLib::IndexValueVector<GlobalIndexType>>
        _dirichlet_bcs;
    std::vector<std::unique_ptr<BoundaryCondition>> _boundary_conditions;
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        _parameters;
};

}  // namespace ProcessLib
