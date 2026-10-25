// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/ProcessVariable.h"
#include "SourceTerm.h"

namespace ProcessLib
{
class SourceTermCollection final
{
public:
    explicit SourceTermCollection(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters)
        : _parameters(parameters)
    {
    }

    void integrate(const double t, GlobalVector const& x, GlobalVector& b,
                   GlobalMatrix* jac) const;

    void updateActiveElements(std::vector<std::size_t> const& ids,bool masked);
    bool hasHistory() const;
    bool supportsCheckpoint() const;
    bool supportsSafety() const;
    void setSafetyState(bool active,double factor,double physical_time);
    std::vector<std::vector<double>> checkpointState() const;
    void restoreCheckpointState(std::vector<std::vector<double>> const&);

    void postAcceptedTimestep(double t, GlobalVector const& x);

    void addSourceTermsForProcessVariables(
        std::vector<std::reference_wrapper<ProcessVariable>> const&
            process_variables,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        unsigned const integration_order,
        const MeshLib::Mesh& bulk_mesh);

private:
    bool safety_active_=false;
    double safety_physical_time_=0.,safety_factor_=1.;
    std::vector<std::unique_ptr<SourceTermBase>> _source_terms;
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        _parameters;
};

}  // namespace ProcessLib
