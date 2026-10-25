// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "SourceTerm.h"
#include "ProcessLib/StructuralElements/StructuralKernels.h"

namespace ProcessLib
{
template <int GlobalDim>
class AnchorTerm final : public SourceTerm
{
public:
    explicit AnchorTerm(
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
        std::size_t const source_term_mesh_id, MeshLib::Mesh const& st_mesh,
        const int variable_id,
        ParameterLib::Parameter<double> const& parameter,
        bool plastic = false, double tension_limit = 0.,
        double compression_limit = 0., bool tension_only = false, double initial_force=0., double activation_time=0., double prestress_end_time=-1.);

    void integrate(const double t, GlobalVector const& x, GlobalVector& b,
                   GlobalMatrix* jac) const override;

    bool hasHistory() const override { return plastic_ || initial_force_!=0. || activation_time_>0. || prestress_end_time_>=0.; }

    bool supportsSafety() const override { return true; }
    bool supportsCheckpoint() const override { return true; }
    std::vector<double> checkpointState() const override
    {
        if(!hasHistory()) { return {}; }
        std::vector<double> result{2.,last_observed_time_,static_cast<double>(accepted_states_.size())};
        result.insert(result.end(),reference_strains_.begin(),reference_strains_.end());
        for(auto const& state:accepted_states_) { result.push_back(state.plastic_strain); }
        return result;
    }
    void restoreCheckpointState(std::vector<double> const& state) override
    {
        if(!hasHistory() && state.empty()) { return; }
        auto const n=accepted_states_.size();
        if(state.size()!=3+2*n || state[0]!=2. || state[2]!=n || state[1]<0. ||
            !std::all_of(state.begin(),state.end(),[](double v){return std::isfinite(v);}))
        { OGS_FATAL("Invalid anchor checkpoint state."); }
        last_observed_time_=state[1];
        for(std::size_t i=0;i<n;++i) { reference_strains_[i]=state[3+i];accepted_states_[i].plastic_strain=state[3+n+i]; }
    }

    void postAcceptedTimestep(double t, GlobalVector const& x) override;

private:
    bool inactive(double t) const { return activation_time_>0. && t<=activation_time_+1e-12*std::max(1.,activation_time_); }
    bool forceControlled(double t) const { return prestress_end_time_>=0. && !inactive(t) && t<=prestress_end_time_+1e-12*std::max(1.,prestress_end_time_); }
    double const prestress_end_time_;
    double const initial_force_;
    double const activation_time_;
    double last_observed_time_=0.;
    std::vector<double> reference_strains_;
    bool const plastic_;
    double const tension_limit_;
    double const compression_limit_;
    bool const tension_only_;
    std::vector<StructuralElements::AxialState> accepted_states_;
    std::size_t const source_term_mesh_id_;
    MeshLib::Mesh const& st_mesh_;
    int const variable_id_;
    ParameterLib::Parameter<double> const& parameter_;
};

extern template class AnchorTerm<2>;
extern template class AnchorTerm<3>;
}  // namespace ProcessLib
