// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <Eigen/Core>
#include "Beam2DProcessData.h"
#include "StructuralKernels.h"
namespace ProcessLib::StructuralElements
{
/// Elastic two-node planar frame, q=[ux1,uy1,theta1,ux2,uy2,theta2].
class Beam2DLocalAssembler
{
public:
    Beam2DLocalAssembler(Eigen::Vector2d const& delta,
                         Beam2DProcessData const& data)
    {
        double const L=delta.norm();
        local_=frameStiffness(L,data.EA,data.EI,data.GA);
        double const c=delta[0]/L,s=delta[1]/L;
        Eigen::Matrix3d R;
        R<<c,s,0.,-s,c,0.,0.,0.,1.;
        transform_.setZero();
        transform_.block<3,3>(0,0)=R;
        transform_.block<3,3>(3,3)=R;
        stiffness_=transform_.transpose()*local_*transform_;
    }
    Eigen::Matrix<double,6,6> const& stiffness() const { return stiffness_; }
    Eigen::Vector<double,6> localForces(Eigen::Vector<double,6> const& q) const
    { return local_*(transform_*q); }
private:
    Eigen::Matrix<double,6,6> local_,transform_,stiffness_;
};
}
