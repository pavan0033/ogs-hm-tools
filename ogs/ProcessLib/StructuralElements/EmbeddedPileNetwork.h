// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <Eigen/QR>
#include <optional>
#include <algorithm>
#include <vector>
#include "Beam2DLocalAssembler.h"
#include "SectionAndInterfaceKernels.h"
namespace ProcessLib::StructuralElements
{
struct PileProperties
{
    double EA,EI,shaft_stiffness,shaft_limit,lateral_stiffness,lateral_limit;
    double base_stiffness,base_limit;
};
struct PileControls
{
    Eigen::VectorXd external_force;
    std::optional<double> axial_displacement,lateral_displacement;
    bool bond_head_to_soil=false;
    std::optional<double> prestress_force;
};
struct PileState
{
    Eigen::VectorXd q,reference;
    std::vector<InterfaceState> shaft,lateral;
    InterfaceState base;
    double prestress_strain=0.;
};
struct PileResponse
{
    PileState state;
    Eigen::VectorXd interface_force,shaft_force,lateral_force;
    double base_force;
    Eigen::VectorXd head_reaction;
    Eigen::MatrixXd soil_tangent;
};
/// Independent elastic beam translations/rotations with distributed elastic-
/// plastic shaft/lateral slip and unilateral end bearing. All trials are pure.
class EmbeddedPileNetwork
{
public:
    EmbeddedPileNetwork(std::vector<double> const& s,
                         std::vector<std::array<std::size_t,2>> const& members,
                         std::size_t head,std::size_t tip,PileProperties properties,
                         std::vector<double> interaction_factors={},
                         std::vector<std::array<double,2>> member_stiffness={})
        :p_(properties),head_(head),tip_(tip),weights_(Eigen::VectorXd::Zero(s.size())),
         K_(Eigen::MatrixXd::Zero(3*s.size(),3*s.size()))
    {
        if(s.size()<2 || head>=s.size() || tip>=s.size() || head==tip ||
           !std::isfinite(p_.shaft_stiffness) || p_.shaft_stiffness<0. ||
           !std::isfinite(p_.shaft_limit) || p_.shaft_limit<0. ||
           !std::isfinite(p_.lateral_stiffness) || p_.lateral_stiffness<0. ||
           std::isnan(p_.lateral_limit) || p_.lateral_limit<0. ||
           !std::isfinite(p_.base_stiffness) || p_.base_stiffness<0. ||
           !std::isfinite(p_.base_limit) || p_.base_limit<0.)
        { throw std::invalid_argument("Invalid embedded pile properties"); }
        if((!interaction_factors.empty() && interaction_factors.size()!=members.size()) ||
           (!member_stiffness.empty() && member_stiffness.size()!=members.size()))
        { throw std::invalid_argument("Invalid pile member profile sizes"); }
        prestress_vector_=Eigen::VectorXd::Zero(3*s.size());
        std::vector<bool> used(s.size(),false);std::size_t member_id=0;
        for(auto const& e:members)
        {
            double const transfer=interaction_factors.empty()?1.:interaction_factors[member_id];
            if(!std::isfinite(transfer) || transfer<0. || transfer>1.)
            { throw std::invalid_argument("Pile interaction factors must be in [0,1]"); }
            auto const stiffness=member_stiffness.empty()?std::array<double,2>{p_.EA,p_.EI}:member_stiffness[member_id];
            ++member_id;used.at(e[0])=true;used.at(e[1])=true;
            double const ds=s.at(e[1])-s.at(e[0]);
            if(transfer==0.)
            {
                double const sign=ds>0.?1.:-1.;
                prestress_vector_[3*e[0]]-=sign*stiffness[0];prestress_vector_[3*e[1]]+=sign*stiffness[0];
                prestress_denominator_+=stiffness[0]*std::abs(ds);prestress_length_+=std::abs(ds);
            }
            if(!std::isfinite(ds) || ds==0.) { throw std::invalid_argument("Embedded member length must be finite and nonzero"); }
            Eigen::Matrix<double,6,6> member_matrix;
            if(stiffness[1]==0.)
            {
                if(!std::isfinite(stiffness[0]) || stiffness[0]<=0.) { throw std::invalid_argument("Embedded rod EA must be positive finite"); }
                member_matrix.setZero();double const k=stiffness[0]/std::abs(ds);
                member_matrix(0,0)=k;member_matrix(3,3)=k;member_matrix(0,3)=-k;member_matrix(3,0)=-k;
            }
            else
            {
                Beam2DLocalAssembler const member(Eigen::Vector2d(ds,0.),{stiffness[0],stiffness[1],std::numeric_limits<double>::infinity()});
                member_matrix=member.stiffness();
            }
            weights_[e[0]]+=transfer*std::abs(ds)/2.;weights_[e[1]]+=transfer*std::abs(ds)/2.;
            for(int a=0;a<2;++a)for(int b=0;b<2;++b)
            { K_.block<3,3>(3*e[a],3*e[b])+=member_matrix.block<3,3>(3*a,3*b); }
        }
        if(std::find(used.begin(),used.end(),false)!=used.end() || weights_.sum()<=0.) { throw std::invalid_argument("Unused pile nodes or no interacting length"); }
    }
    PileState initialState() const
    {
        return {Eigen::VectorXd::Zero(K_.rows()),Eigen::VectorXd::Zero(K_.rows()),
                std::vector<InterfaceState>(weights_.size()),std::vector<InterfaceState>(weights_.size()),{}};
    }
    Eigen::VectorXd const& weights() const { return weights_; }
    double prestressLength() const { return prestress_length_; }
    PileResponse trial(Eigen::VectorXd const& soil,PileState const& accepted,
                        PileControls const& controls,double resistance_factor=1.) const
    {
        if(!std::isfinite(resistance_factor) || resistance_factor<1. || soil.size()!=K_.rows() || controls.external_force.size()!=K_.rows() ||
            accepted.q.size()!=K_.rows() || accepted.reference.size()!=K_.rows() ||
            accepted.shaft.size()!=static_cast<std::size_t>(weights_.size()) ||
            accepted.lateral.size()!=accepted.shaft.size() || !soil.allFinite() ||
            !controls.external_force.allFinite() || !accepted.q.allFinite() || !accepted.reference.allFinite() ||
            (controls.axial_displacement && !std::isfinite(*controls.axial_displacement)) ||
            (controls.lateral_displacement && !std::isfinite(*controls.lateral_displacement)))
        { throw std::invalid_argument("Invalid pile trial vectors or controls"); }
        if(!std::isfinite(accepted.prestress_strain) ||
           (controls.prestress_force && (!std::isfinite(*controls.prestress_force) || prestress_denominator_<=0.)))
        { throw std::invalid_argument("Embedded prestressing requires finite force and a free interaction length"); }
        if(controls.bond_head_to_soil && (controls.axial_displacement || controls.lateral_displacement))
        { throw std::invalid_argument("Bonded pile head cannot have independent displacement controls"); }
        auto q=accepted.q;std::vector<Eigen::Index> free;
        std::vector<bool> constrained(q.size(),false);
        if(controls.axial_displacement) { constrained[3*head_]=true;q[3*head_]=accepted.reference[3*head_]+*controls.axial_displacement; }
        if(controls.lateral_displacement) { constrained[3*head_+1]=true;q[3*head_+1]=accepted.reference[3*head_+1]+*controls.lateral_displacement; }
        if(controls.bond_head_to_soil)
        { for(int c=0;c<2;++c) { constrained[3*head_+c]=true;q[3*head_+c]=soil[3*head_+c]; } }
        for(Eigen::Index i=0;i<q.size();++i)if(!constrained[i]) { free.push_back(i); }
        Eigen::MatrixXd elastic=K_;
        for(Eigen::Index i=0;i<weights_.size();++i)
        { elastic(3*i,3*i)+=p_.shaft_stiffness*weights_[i];elastic(3*i+1,3*i+1)+=p_.lateral_stiffness*weights_[i]; }
        elastic(3*tip_,3*tip_)+=p_.base_stiffness;
        auto submatrix=[&](Eigen::MatrixXd const& matrix)
        {
            Eigen::MatrixXd out(free.size(),free.size());
            for(std::size_t i=0;i<free.size();++i)for(std::size_t j=0;j<free.size();++j)out(i,j)=matrix(free[i],free[j]);
            return out;
        };
        for(int iteration=0;iteration<100;++iteration)
        {
            auto r=assemble(q,soil,accepted,controls.external_force,resistance_factor,controls.prestress_force);
            Eigen::VectorXd g(free.size());for(std::size_t i=0;i<free.size();++i) { g[i]=r.residual[free[i]]; }
            Eigen::MatrixXd const ff=submatrix(r.tangent);
            if(g.lpNorm<Eigen::Infinity>()<1e-8)
            {
                Eigen::MatrixXd D=Eigen::MatrixXd::Zero(q.size(),q.size());D.diagonal()=r.diagonal;
                Eigen::MatrixXd df(free.size(),q.size());
                for(std::size_t i=0;i<free.size();++i) { df.row(i)=D.row(free[i]); }
                if(controls.bond_head_to_soil)
                {
                    Eigen::MatrixXd H=Eigen::MatrixXd::Zero(q.size(),q.size());
                    H(3*head_,3*head_)=1.;H(3*head_+1,3*head_+1)=1.;
                    Eigen::MatrixXd const prescribed=r.tangent*H;
                    for(std::size_t i=0;i<free.size();++i) { df.row(i)-=prescribed.row(free[i]); }
                    Eigen::MatrixXd const free_recovery=ff.completeOrthogonalDecomposition().solve(df);
                    Eigen::MatrixXd Q=H;
                    for(std::size_t i=0;i<free.size();++i) { Q.row(free[i])=free_recovery.row(i); }
                    r.interface_force-=H.transpose()*r.residual;
                    Eigen::MatrixXd tangent=D-D*Q+H.transpose()*(r.tangent*Q-D);
                    return {PileState{q,accepted.reference,std::move(r.shaft),std::move(r.lateral),r.base,controls.prestress_force?(prestress_vector_.dot(q-accepted.reference)-*controls.prestress_force*prestress_length_)/prestress_denominator_:accepted.prestress_strain},
                            std::move(r.interface_force),std::move(r.shaft_force),std::move(r.lateral_force),r.base_force,
                            std::move(r.residual),std::move(tangent)};
                }
                Eigen::MatrixXd const recovery=ff.completeOrthogonalDecomposition().solve(df);
                return {PileState{q,accepted.reference,std::move(r.shaft),std::move(r.lateral),r.base,controls.prestress_force?(prestress_vector_.dot(q-accepted.reference)-*controls.prestress_force*prestress_length_)/prestress_denominator_:accepted.prestress_strain},
                        std::move(r.interface_force),std::move(r.shaft_force),std::move(r.lateral_force),r.base_force,
                        std::move(r.residual),D-df.transpose()*recovery};
            }
            Eigen::VectorXd step=ff.completeOrthogonalDecomposition().solve(-g);
            if(!step.allFinite() || g.dot(step)>=-1e-12*g.squaredNorm())
            { step=submatrix(elastic).completeOrthogonalDecomposition().solve(-g); }
            if(!step.allFinite()) { break; }
            bool improved=false;
            for(double scale=1.;scale>=1e-9;scale*=.5)
            {
                auto candidate=q;
                for(std::size_t i=0;i<free.size();++i) { candidate[free[i]]+=scale*step[i]; }
                auto const next=assemble(candidate,soil,accepted,controls.external_force,resistance_factor,controls.prestress_force);
                // Near equilibrium the absolute stored energies can round to
                // the same value. Accept a decreasing free-force merit as an
                // alternative, without changing the equilibrium tolerance.
                Eigen::VectorXd next_g(free.size());
                for (std::size_t i=0;i<free.size();++i)
                { next_g[i]=next.residual[free[i]]; }
                bool const force_decrease =
                    next_g.norm() <= (1.-1e-4*scale)*g.norm();
                bool const energy_decrease = next.energy < r.energy &&
                    next.energy<=r.energy+1e-4*scale*g.dot(step);
                if(force_decrease || energy_decrease)
                { q=std::move(candidate);improved=true;break; }
            }
            if(!improved) { break; }
        }
        throw std::runtime_error("Embedded pile internal equilibrium failed; reduce increment or check capacity.");
    }
private:
    struct Assembly
    {
        Eigen::VectorXd residual,interface_force,diagonal,shaft_force,lateral_force;
        Eigen::MatrixXd tangent;
        double energy=0.,base_force=0.;
        std::vector<InterfaceState> shaft,lateral;
        InterfaceState base;
    };
    Assembly assemble(Eigen::VectorXd const& q,Eigen::VectorXd const& soil,
                       PileState const& accepted,Eigen::VectorXd const& external,double resistance_factor,std::optional<double> prestress_force) const
    {
        Eigen::VectorXd const relative=q-soil;
        Eigen::VectorXd const dq=q-accepted.reference;
        Eigen::MatrixXd effective=K_;Eigen::VectorXd offset=-prestress_vector_*accepted.prestress_strain;
        if(prestress_force)
        {
            effective-=prestress_vector_*prestress_vector_.transpose()/prestress_denominator_;
            offset=prestress_vector_*(*prestress_force*prestress_length_/prestress_denominator_);
        }
        Assembly out;// Stiff member forces and prestress nearly cancel at equilibrium.
        // Accumulate this small local product in extended precision so the
        // unchanged absolute residual criterion is not limited by cancellation.
        out.residual =
            (effective.cast<long double>() * dq.cast<long double>() +
             offset.cast<long double>() - external.cast<long double>())
                .cast<double>();out.interface_force=Eigen::VectorXd::Zero(q.size());out.diagonal=out.interface_force;
        out.shaft_force=Eigen::VectorXd::Zero(weights_.size());out.lateral_force=out.shaft_force;
        out.energy=.5*dq.dot(effective*dq)+offset.dot(dq)-external.dot(q);
        out.shaft=accepted.shaft;out.lateral=accepted.lateral;out.base=accepted.base;
        auto add=[&](Eigen::Index dof,double k,double limit,InterfaceState old,bool one_sided,InterfaceState& state)
        {
            if(k==0. || limit==0.) { return 0.; }
            InterfaceResponse r;
            if(std::isinf(limit)) { r={k*relative[dof],k,old}; }
            else { r=slipTrial(relative[dof],old,k,limit,one_sided); }
            out.interface_force[dof]+=r.resistance;out.diagonal[dof]+=r.tangent;state=r.trial_state;
            out.energy+=.5*r.resistance*r.resistance/k;
            if(std::isfinite(limit)) { out.energy+=limit*(one_sided?std::max(0.,state.plastic_slip-old.plastic_slip):std::abs(state.plastic_slip-old.plastic_slip)); }
            return r.resistance;
        };
        for(Eigen::Index i=0;i<weights_.size();++i)
        {
            out.shaft_force[i]=add(3*i,p_.shaft_stiffness*weights_[i],p_.shaft_limit*weights_[i]/resistance_factor,accepted.shaft[i],false,out.shaft[i]);
            out.lateral_force[i]=add(3*i+1,p_.lateral_stiffness*weights_[i],p_.lateral_limit*weights_[i]/resistance_factor,accepted.lateral[i],false,out.lateral[i]);
        }
        out.base_force=add(3*tip_,p_.base_stiffness,p_.base_limit/resistance_factor,accepted.base,true,out.base);
        out.residual+=out.interface_force;out.tangent=effective;out.tangent.diagonal()+=out.diagonal;
        return out;
    }
    PileProperties p_;
    std::size_t head_,tip_;
    Eigen::VectorXd weights_;
    Eigen::MatrixXd K_;
    Eigen::VectorXd prestress_vector_;
    double prestress_denominator_=0.,prestress_length_=0.;
};
}
