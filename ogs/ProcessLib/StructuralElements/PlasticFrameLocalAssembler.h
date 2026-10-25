// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <array>
#include "SectionAndInterfaceKernels.h"

namespace ProcessLib::StructuralElements
{
struct PlasticFrameResponse
{
    double incremental_potential=0.;
    Eigen::Vector<double,6> force=Eigen::Vector<double,6>::Zero();
    Eigen::Matrix<double,6,6> tangent=Eigen::Matrix<double,6,6>::Zero();
    std::array<SectionState,2> states;
    std::array<Eigen::Vector2d,2> section_forces;
};
/// Two-point displacement-based Euler beam; pure section trials from accepted
/// plastic strain/curvature. No geometric nonlinearity or shear deformation.
class PlasticFrameLocalAssembler
{
public:
    PlasticFrameLocalAssembler(Eigen::Vector2d const& delta,double EA,double EI,
                                double Np,double Mp)
        : L_(delta.norm()), stiffness_(EA,EI), capacities_(Np,Mp)
    {
        if(L_<=0. || !std::isfinite(L_))
        { throw std::invalid_argument("Invalid plastic beam length"); }
        double const c=delta[0]/L_,s=delta[1]/L_;
        Eigen::Matrix3d R;
        R<<c,s,0.,-s,c,0.,0.,0.,1.;
        transform_.setZero();transform_.block<3,3>(0,0)=R;transform_.block<3,3>(3,3)=R;
        // Validate all section properties before the first global iteration.
        sectionTrial(Eigen::Vector2d::Zero(),SectionState{},stiffness_,capacities_);
    }
    PlasticFrameResponse trial(Eigen::Vector<double,6> const& q,
                                std::array<SectionState,2> const& accepted) const
    {
        PlasticFrameResponse out;
        Eigen::Vector<double,6> const local=transform_*q;
        for(int i=0;i<2;++i)
        {
            double const xi=(1.+(i==0?-1.:1.)/std::sqrt(3.))/2.;
            Eigen::Matrix<double,2,6> B=Eigen::Matrix<double,2,6>::Zero();
            B(0,0)=-1./L_;B(0,3)=1./L_;
            B(1,1)=(-6.+12.*xi)/(L_*L_);B(1,2)=(-4.+6.*xi)/L_;
            B(1,4)=(6.-12.*xi)/(L_*L_);B(1,5)=(-2.+6.*xi)/L_;
            auto const response=sectionTrial(B*local,accepted[i],stiffness_,capacities_);
            Eigen::Vector2d const increment=response.trial_state.plastic_deformation-accepted[i].plastic_deformation;
            out.incremental_potential+=(.5*(response.forces.array().square()/stiffness_.array()).sum()+
                capacities_.cwiseProduct(increment).cwiseAbs().maxCoeff())*L_/2.;
            out.force.noalias()+=B.transpose()*response.forces*L_/2.;
            out.tangent.noalias()+=B.transpose()*response.tangent*B*L_/2.;
            out.states[i]=response.trial_state;out.section_forces[i]=response.forces;
        }
        out.force=transform_.transpose()*out.force;
        out.tangent=transform_.transpose()*out.tangent*transform_;
        return out;
    }
    Eigen::Vector<double,6> localForces(Eigen::Vector<double,6> const& global) const
    { return transform_*global; }
private:
    double L_;
    Eigen::Vector2d stiffness_,capacities_;
    Eigen::Matrix<double,6,6> transform_;
};
}
