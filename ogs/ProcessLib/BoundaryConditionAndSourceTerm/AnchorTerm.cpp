// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "AnchorTerm.h"

#include <Eigen/Core>
#include <cassert>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
template <int GlobalDim>
AnchorTerm<GlobalDim>::AnchorTerm(
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
    std::size_t const source_term_mesh_id,
    MeshLib::Mesh const& st_mesh,
    const int variable_id,
    ParameterLib::Parameter<double> const& parameter,
    bool const plastic, double const tension_limit,
    double const compression_limit, bool const tension_only,double const initial_force,double const activation_time,double const prestress_end_time)
    : SourceTerm(std::move(source_term_dof_table)),
      prestress_end_time_(prestress_end_time),initial_force_(initial_force),activation_time_(activation_time),
      reference_strains_(st_mesh.getNumberOfElements(),0.),
      plastic_(plastic),
      tension_limit_(tension_limit),
      compression_limit_(compression_limit),
      tension_only_(tension_only),
      accepted_states_(st_mesh.getNumberOfElements()),
      source_term_mesh_id_(source_term_mesh_id),
      st_mesh_(st_mesh),
      variable_id_(variable_id),
      parameter_(parameter)
{
    DBUG("Create AnchorTerm.");
    if(plastic_ && initial_force_!=0.)
    {
        for(auto const* e:st_mesh_.getElements())
        {
            ParameterLib::SpatialPosition pos;pos.setElementID(e->getID());
            double const EA=parameter_(0.,pos)[0];
            if(!std::isfinite(EA) || EA<=0.) { OGS_FATAL("Anchor EA must be positive finite."); }
            accepted_states_[e->getID()].plastic_strain=-initial_force_/EA;
        }
    }
}

template <int GlobalDim>
void AnchorTerm<GlobalDim>::integrate(const double t, GlobalVector const& x,
                                      GlobalVector& b, GlobalMatrix* jac) const
{
    DBUG("Assemble AnchorTerm.");
    if(inactive(t)) { return; }
    if(activation_time_>0. && last_observed_time_<activation_time_-1e-12*std::max(1.,activation_time_))
    { OGS_FATAL("Anchor activation_time must coincide with an accepted time boundary."); }

    if(prestress_end_time_>=0. && t>prestress_end_time_+1e-12*std::max(1.,prestress_end_time_) &&
       last_observed_time_<prestress_end_time_-1e-12*std::max(1.,prestress_end_time_))
    { OGS_FATAL("Anchor prestress_end_time must coincide with an accepted time boundary."); }

    using GlobalDimVector = Eigen::Vector<double, GlobalDim>;
    using GlobalDimMatrix = Eigen::Matrix<double, GlobalDim, GlobalDim>;

    for (MeshLib::Element const* const element : st_mesh_.getElements())
    {
        auto const element_id = element->getID();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(element_id);

        std::vector<GlobalIndexType> const global_indices =
            NumLib::getIndices(element_id, *_source_term_dof_table);
        assert(global_indices.size() == 2 * GlobalDim);

        Eigen::Vector<double, 2 * GlobalDim> const local_x =
            MathLib::toVector(x.get(global_indices));

        Eigen::Vector<double, 2 * GlobalDim> local_rhs =
            Eigen::Vector<double, 2 * GlobalDim>::Zero();
        Eigen::Matrix<double, 2 * GlobalDim, 2 * GlobalDim> local_Jac =
            Eigen::Matrix<double, 2 * GlobalDim, 2 * GlobalDim>::Zero();

        // The local indices are, due to the nature of the DOF table, all even
        // for the first node and odd for the second node.
        auto node_local_indices = [](int const i)
        { return Eigen::seqN(i, Eigen::fix<GlobalDim>, Eigen::fix<2>); };

        auto node_coords = [element](int const i)
        { return element->getNode(i)->asEigenVector3d(); };
        GlobalDimVector const l_original =
            (node_coords(1) - node_coords(0)).template head<GlobalDim>();
        double const l_original_norm = l_original.norm();

        // Displacement in the two nodes.
        auto u = [&local_x, &node_local_indices](int const i)
        { return local_x(node_local_indices(i)); };
        GlobalDimVector const l = l_original + u(1) - u(0);

        double const K = parameter_(t, pos)[0];
        GlobalDimVector f;
        GlobalDimMatrix Df;
        if (forceControlled(t))
        {
            f=l_original/l_original_norm*initial_force_;
            Df.setZero();
        }
        else if (plastic_)
        {
            GlobalDimVector const direction = l_original / l_original_norm;
            double const strain = direction.dot(u(1)-u(0)) / l_original_norm-reference_strains_[element_id];
            auto const response = StructuralElements::axialTrial(
                strain, accepted_states_.at(element_id), K, tension_limit_,
                compression_limit_, tension_only_);
            f = direction * response.force;
            Df = response.tangent / l_original_norm * direction * direction.transpose();
        }
        else
        {
            f = l_original / l_original_norm *
                (K*((l.norm() - l_original_norm) / l_original_norm-reference_strains_[element_id])+initial_force_);
            Df = l_original / l_original_norm * K *
                 l.transpose() / l.norm() / l_original_norm;
        }

        // Signs for the two nodes alternate.
        constexpr auto even_odd_sign = [](int const n)
        { return (n % 2 == 0) ? 1.0 : -1.0; };

        // There are two blocks in the rhs vector for the two nodes of the
        // element and four blocks in the local Jacobian formed by the four
        // combinations of the two nodes.
        for (int i = 0; i < 2; ++i)
        {
            local_rhs(node_local_indices(i)).noalias() += even_odd_sign(i) * f;

            for (int j = 0; j < 2; ++j)
            {
                local_Jac(node_local_indices(i), node_local_indices(j))
                    .noalias() += even_odd_sign(i) * even_odd_sign(j) * Df;
            }
        }

        b.add(global_indices, local_rhs);
        if (jac)
        {
            jac->add({global_indices, global_indices}, local_Jac);
        }
    }
}

template <int GlobalDim>
void AnchorTerm<GlobalDim>::postAcceptedTimestep(double const t, GlobalVector const& x)
{
    if (!hasHistory()) { return; }
    for (auto const* element : st_mesh_.getElements())
    {
        auto const id = element->getID();
        auto const indices = NumLib::getIndices(id, *_source_term_dof_table);
        auto const values = x.get(indices);
        Eigen::Vector<double, GlobalDim> const delta =
            (element->getNode(1)->asEigenVector3d() -
             element->getNode(0)->asEigenVector3d()).template head<GlobalDim>();
        Eigen::Vector<double, GlobalDim> du;
        for (int c=0; c<GlobalDim; ++c) { du[c]=values[2*c+1]-values[2*c]; }
        ParameterLib::SpatialPosition pos;
        pos.setElementID(id);
        double const observed=plastic_?delta.dot(du)/delta.squaredNorm():((delta+du).norm()-delta.norm())/delta.norm();
        if(inactive(t)) { reference_strains_[id]=observed;continue; }
        if(forceControlled(t))
        {
            if(plastic_) { accepted_states_[id].plastic_strain=observed-reference_strains_[id]-initial_force_/parameter_(t,pos)[0]; }
            else { reference_strains_[id]=observed; }
            INFO("Accepted structural element {}: t={}, force={}, plastic_strain={}",id,t,initial_force_,accepted_states_[id].plastic_strain);
            continue;
        }
        if(!plastic_)
        {
            double const force=parameter_(t,pos)[0]*(observed-reference_strains_[id])+initial_force_;
            INFO("Accepted structural element {}: t={}, force={}, plastic_strain={}",id,t,force,0.);
            continue;
        }
        auto const response = StructuralElements::axialTrial(
            observed-reference_strains_[id], accepted_states_.at(id),
            parameter_(t,pos)[0], tension_limit_, compression_limit_, tension_only_);
        accepted_states_.at(id) = response.trial_state;
        INFO("Accepted structural element {}: t={}, force={}, plastic_strain={}",
             id, t, response.force, response.trial_state.plastic_strain);
    }
    last_observed_time_=t;
}

template class AnchorTerm<2>;
template class AnchorTerm<3>;

}  // namespace ProcessLib
