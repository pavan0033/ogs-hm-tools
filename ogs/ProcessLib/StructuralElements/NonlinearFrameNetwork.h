// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <Eigen/QR>
#include <vector>
#include "PlasticFrameLocalAssembler.h"
#include "Beam2DLocalAssembler.h"

namespace ProcessLib::StructuralElements
{
struct FrameNetworkState
{
    Eigen::VectorXd displacement;
    Eigen::VectorXd rotations;
    std::vector<std::array<SectionState,2>> sections;
};
struct FrameNetworkResponse
{
    Eigen::VectorXd force;
    Eigen::MatrixXd tangent;
    FrameNetworkState state;
    std::vector<Eigen::Vector<double,6>> local_end_forces;
};
/// Pure nonlinear condensation of shared unloaded rotations. Accepted history is
/// supplied explicitly and is never modified by a trial or failed local solve.
class NonlinearFrameNetwork
{
public:
    NonlinearFrameNetwork(std::vector<Eigen::Vector2d> const& points,
                          std::vector<std::array<std::size_t,2>> members,
                          double EA,double EI,double Np,double Mp, Eigen::VectorXd external_load=Eigen::VectorXd{})
        : n_(points.size()), members_(std::move(members))
    {
        load_=external_load.size()==0?Eigen::VectorXd::Zero(3*n_).eval():std::move(external_load);
        if(load_.size()!=3*n_ || !load_.allFinite()) { throw std::invalid_argument("Invalid plate distributed load vector"); }
        elastic_rr_=Eigen::MatrixXd::Zero(n_,n_);
        for(auto const& edge:members_)
        {
            Eigen::Vector2d const delta=points.at(edge[1])-points.at(edge[0]);
            assemblers_.emplace_back(delta,EA,EI,Np,Mp);
            Beam2DLocalAssembler const elastic(delta,{EA,EI,std::numeric_limits<double>::infinity()});
            for(int a=0;a<2;++a) for(int b=0;b<2;++b)
            { elastic_rr_(edge[a],edge[b])+=elastic.stiffness()(3*a+2,3*b+2); }
        }
    }
    FrameNetworkState initialState() const
    {
        return {Eigen::VectorXd::Zero(2*n_),Eigen::VectorXd::Zero(n_),
                std::vector<std::array<SectionState,2>>(members_.size())};
    }
    FrameNetworkResponse trial(Eigen::VectorXd const& u,
                                FrameNetworkState const& accepted,
                                Eigen::VectorXd theta) const
    {
        for(int iteration=0;iteration<80;++iteration)
        {
            auto response=assemble(u,theta,accepted);
            
            if(response.r.lpNorm<Eigen::Infinity>()<1e-8)
            {
                auto const recovery=response.rr.completeOrthogonalDecomposition().solve(response.tr.transpose()).eval();
                return {std::move(response.f),response.tt-response.tr*recovery,
                        FrameNetworkState{u,theta,std::move(response.sections)},
                        std::move(response.local_end_forces)};
            }
            Eigen::VectorXd delta=response.rr.completeOrthogonalDecomposition().solve(-response.r);
            if(!delta.allFinite() || response.r.dot(delta)>=-1e-12*response.r.squaredNorm())
            { delta=elastic_rr_.ldlt().solve(-response.r); }
            if(!delta.allFinite()) { break; }
            bool improved=false;
            for(double scale=1.;scale>=1e-8;scale*=.5)
            {
                Eigen::VectorXd candidate=theta+scale*delta;
                auto const next=assemble(u,candidate,accepted);
                if(next.energy<=response.energy+1e-4*scale*response.r.dot(delta))
                {
                    theta=std::move(candidate);improved=true;break;
                }
            }
            if(!improved) { break; }
        }
        throw std::runtime_error("Nonlinear plate rotation equilibrium failed; reduce global increment.");
    }
private:
    struct Assembly
    {
        double energy=0.;
        Eigen::VectorXd f,r;
        Eigen::MatrixXd tt,tr,rr;
        std::vector<std::array<SectionState,2>> sections;
        std::vector<Eigen::Vector<double,6>> local_end_forces;
    };
    Assembly assemble(Eigen::VectorXd const& u,Eigen::VectorXd const& theta,
                      FrameNetworkState const& accepted) const
    {
        Assembly out{0.,Eigen::VectorXd::Zero(2*n_),Eigen::VectorXd::Zero(n_),
                     Eigen::MatrixXd::Zero(2*n_,2*n_),Eigen::MatrixXd::Zero(2*n_,n_),
                     Eigen::MatrixXd::Zero(n_,n_),{}, {}};
        for(std::size_t e=0;e<members_.size();++e)
        {
            Eigen::Vector<double,6> q;
            for(int a=0;a<2;++a)
            {
                auto const i=members_[e][a];q.segment<2>(3*a)=u.segment<2>(2*i);q[3*a+2]=theta[i];
            }
            auto const response=assemblers_[e].trial(q,accepted.sections.at(e));
            out.energy+=response.incremental_potential;
            out.sections.push_back(response.states);
            out.local_end_forces.push_back(assemblers_[e].localForces(response.force));
            for(int a=0;a<2;++a)
            {
                auto const i=members_[e][a];
                out.f.segment<2>(2*i)+=response.force.segment<2>(3*a);
                out.r[i]+=response.force[3*a+2];
                for(int b=0;b<2;++b)
                {
                    auto const j=members_[e][b];
                    out.tt.block<2,2>(2*i,2*j)+=response.tangent.block<2,2>(3*a,3*b);
                    out.tr.block<2,1>(2*i,j)+=response.tangent.block<2,1>(3*a,3*b+2);
                    out.rr(i,j)+=response.tangent(3*a+2,3*b+2);
                }
            }
        }
        for(std::size_t i=0;i<n_;++i)
        {
            out.f.segment<2>(2*i)-=load_.segment<2>(3*i);
            out.r[i]-=load_[3*i+2];
            out.energy-=load_.segment<2>(3*i).dot(u.segment<2>(2*i))+load_[3*i+2]*theta[i];
        }
        return out;
    }
    Eigen::VectorXd load_;
    Eigen::MatrixXd elastic_rr_;
    std::size_t n_;
    std::vector<std::array<std::size_t,2>> members_;
    std::vector<PlasticFrameLocalAssembler> assemblers_;
};
}
