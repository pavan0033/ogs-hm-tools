// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ProcessLib::StructuralElements
{
struct AxialState
{
    double plastic_strain = 0.;
};
struct AxialResponse
{
    double force;
    double tangent;
    AxialState trial_state;
};
/// Pure trial evaluation: only the caller may commit trial_state after acceptance.
inline AxialResponse axialTrial(double strain, AxialState const& accepted,
                                double EA, double tension_limit,
                                double compression_limit, bool tension_only)
{
    if (!std::isfinite(strain) || !std::isfinite(accepted.plastic_strain) ||
        !std::isfinite(EA) || EA<=0 || !std::isfinite(tension_limit) || tension_limit<0 ||
        !std::isfinite(compression_limit) || compression_limit<0)
    { throw std::invalid_argument("Invalid axial material input"); }
    double const elastic_strain = strain-accepted.plastic_strain;
    if (tension_only && elastic_strain<=0.) { return {0.,0.,accepted}; }
    double const trial = EA*elastic_strain;
    double const force = std::clamp(trial,-compression_limit,tension_limit);
    if (force==trial) { return {force,EA,accepted}; }
    return {force,0.,AxialState{strain-force/EA}};
}
/// Local 2-node Timoshenko frame: [u1,v1,theta1,u2,v2,theta2].
/// GA includes the shear correction factor. GA=infinity gives Euler-Bernoulli.
inline Eigen::Matrix<double,6,6> frameStiffness(double L,double EA,double EI,double GA)
{
    if (!std::isfinite(L) || L<=0 || !std::isfinite(EA) || EA<=0 ||
        !std::isfinite(EI) || EI<=0 || std::isnan(GA) || GA<=0)
    { throw std::invalid_argument("Invalid frame properties"); }
    double const phi=12.*EI/(GA*L*L);
    double const b=EI/(L*L*L*(1.+phi));
    Eigen::Matrix<double,6,6> K=Eigen::Matrix<double,6,6>::Zero();
    K(0,0)=K(3,3)=EA/L;K(0,3)=K(3,0)=-EA/L;
    Eigen::Matrix4d B;
    B << 12.,6.*L,-12.,6.*L,
         6.*L,(4.+phi)*L*L,-6.*L,(2.-phi)*L*L,
         -12.,-6.*L,12.,-6.*L,
         6.*L,(2.-phi)*L*L,-6.*L,(4.+phi)*L*L;
    int const ids[4]={1,2,4,5};
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)K(ids[i],ids[j])=b*B(i,j);
    return K;
}
}
