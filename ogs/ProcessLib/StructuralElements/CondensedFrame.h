// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <stdexcept>

namespace ProcessLib::StructuralElements
{
/// Exact network-level elimination of unloaded internal rotations.
/// Translation order is [ux0,uy0,...], rotation order is [theta0,...].
struct CondensedFrame
{
    Eigen::MatrixXd stiffness;
    Eigen::MatrixXd rotation_from_translation;
};
inline CondensedFrame condenseFrameRotations(Eigen::MatrixXd const& full)
{
    if (full.rows() != full.cols() || full.rows() == 0 || full.rows() % 3)
    {
        throw std::invalid_argument("Invalid full frame matrix");
    }
    Eigen::Index const n = full.rows() / 3;
    Eigen::MatrixXd tt(2*n,2*n), tr(2*n,n), rr(n,n);
    for (Eigen::Index i=0; i<n; ++i)
    {
        for (Eigen::Index j=0; j<n; ++j)
        {
            tt.block<2,2>(2*i,2*j)=full.block<2,2>(3*i,3*j);
            tr.block<2,1>(2*i,j)=full.block<2,1>(3*i,3*j+2);
            rr(i,j)=full(3*i+2,3*j+2);
        }
    }
    Eigen::LDLT<Eigen::MatrixXd> factor(rr);
    if (factor.info()!=Eigen::Success || !factor.isPositive())
    {
        throw std::invalid_argument("Plate rotation block is not positive definite");
    }
    Eigen::MatrixXd recovery=-factor.solve(tr.transpose());
    return {tt+tr*recovery,recovery};
}
}
