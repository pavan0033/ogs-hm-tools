// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "SourceTerm.h"
#include "ProcessLib/StructuralElements/Beam2DLocalAssembler.h"
#include "ProcessLib/StructuralElements/CondensedFrame.h"
#include "ProcessLib/StructuralElements/NonlinearFrameNetwork.h"
#include "NumLib/Exceptions.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <optional>
#include "ProcessLib/StructuralElements/DistributedLoadRamp.h"

namespace ProcessLib
{
/// Bonded, elastic plate network. Internal rotations are shared at network nodes
/// and eliminated exactly, not independently at each element.
class Plate2D final : public SourceTerm
{
public:
    Plate2D(std::unique_ptr<NumLib::LocalToGlobalIndexMap> dofs,
            MeshLib::Mesh const& mesh, double EA, double EI, double GA, double Np=0.,double Mp=0., Eigen::Vector2d const distributed_load=Eigen::Vector2d::Zero(),
            std::optional<StructuralElements::DistributedLoadRamp> load_ramp=std::nullopt)
        : SourceTerm(std::move(dofs)), mesh_(mesh), data_{EA,EI,GA}, load_ramp_(std::move(load_ramp))
    {
#ifdef USE_PETSC
        OGS_FATAL("Plate2D currently requires serial Eigen.");
#endif
        auto const n=mesh_.getNumberOfNodes();
        if (n==0 || n>1000 || !std::isfinite(EA) || EA<=0. ||
            !std::isfinite(EI) || EI<=0. || std::isnan(GA) || GA<=0.)
        { OGS_FATAL("Plate2D needs positive stiffnesses and 1..1000 nodes."); }
        if(load_ramp_ && Np>0.)
        { OGS_FATAL("Time-ramped distributed loads currently require elastic Plate2D."); }
        if(!distributed_load.allFinite()) { OGS_FATAL("Plate2D distributed_load must be finite."); }
        Eigen::MatrixXd full=Eigen::MatrixXd::Zero(3*n,3*n);
        Eigen::VectorXd load=Eigen::VectorXd::Zero(3*n);
        member_loads_.resize(mesh_.getNumberOfElements());
        indices_.resize(2*n);
        std::vector<bool> used(n,false);
        for (auto const* e:mesh_.getElements())
        {
            if (e->getNumberOfNodes()!=2)
            { OGS_FATAL("Plate2D requires two-node line elements."); }
            auto const delta=(e->getNode(1)->asEigenVector3d()-e->getNode(0)->asEigenVector3d()).eval();
            if (delta.head<2>().norm()<=0. || std::abs(delta[2])>1e-12)
            { OGS_FATAL("Plate2D requires nonzero planar members."); }
            double const length=delta.head<2>().norm();
            Eigen::Vector2d const direction=delta.head<2>()/length;
            Eigen::Vector2d const normal(-direction[1],direction[0]);
            double const axial=distributed_load.dot(direction), transverse=distributed_load.dot(normal);
            Eigen::Vector<double,6> member_load;
            member_load<<axial*length/2,transverse*length/2,transverse*length*length/12,
                         axial*length/2,transverse*length/2,-transverse*length*length/12;
            member_loads_[e->getID()]=member_load;
            for(int a=0;a<2;++a)
            {
                auto const i=e->getNode(a)->getID();
                load.segment<2>(3*i)+=distributed_load*length/2;
                load[3*i+2]+=member_load[3*a+2];
            }
            StructuralElements::Beam2DLocalAssembler const local(delta.head<2>(),data_);
            auto const ids=NumLib::getIndices(e->getID(),*_source_term_dof_table);
            if (ids.size()!=4) { OGS_FATAL("Plate2D needs two displacement components."); }
            for (int a=0;a<2;++a)
            {
                auto const i=e->getNode(a)->getID();used.at(i)=true;
                indices_[2*i]=ids[a];indices_[2*i+1]=ids[2+a];
                for (int b=0;b<2;++b)
                {
                    auto const j=e->getNode(b)->getID();
                    full.block<3,3>(3*i,3*j)+=local.stiffness().block<3,3>(3*a,3*b);
                }
            }
        }
        if (std::find(used.begin(),used.end(),false)!=used.end())
        { OGS_FATAL("Plate2D has unused mesh nodes."); }
        frame_=StructuralElements::condenseFrameRotations(full);
        Eigen::MatrixXd rr(n,n),tr(2*n,n);Eigen::VectorXd ft(2*n),fr(n);
        for(std::size_t i=0;i<n;++i)
        {
            ft.segment<2>(2*i)=load.segment<2>(3*i);fr[i]=load[3*i+2];
            for(std::size_t j=0;j<n;++j)
            { rr(i,j)=full(3*i+2,3*j+2);tr.block<2,1>(2*i,j)=full.block<2,1>(3*i,3*j+2); }
        }
        rotation_offset_=rr.ldlt().solve(fr);
        condensed_load_=ft-tr*rotation_offset_;
        if(Np>0.)
        {
            std::vector<Eigen::Vector2d> points;
            std::vector<std::array<std::size_t,2>> members;
            for(auto const* node:mesh_.getNodes()) { points.push_back(node->asEigenVector3d().head<2>()); }
            for(auto const* e:mesh_.getElements()) { members.push_back({e->getNode(0)->getID(),e->getNode(1)->getID()}); }
            nonlinear_=std::make_unique<StructuralElements::NonlinearFrameNetwork>(points,std::move(members),EA,EI,Np,Mp,load);
            accepted_=nonlinear_->initialState();
        }
    }
    bool hasHistory() const override { return nonlinear_!=nullptr; }
    bool supportsSafety() const override { return !load_ramp_; }
    bool supportsCheckpoint() const override { return true; }
    std::vector<double> checkpointState() const override
    {
        if(!nonlinear_) { return {}; }
        std::vector<double> result{1.,static_cast<double>(accepted_.rotations.size()),static_cast<double>(accepted_.sections.size())};
        result.insert(result.end(),accepted_.displacement.data(),accepted_.displacement.data()+accepted_.displacement.size());
        result.insert(result.end(),accepted_.rotations.data(),accepted_.rotations.data()+accepted_.rotations.size());
        for(auto const& element:accepted_.sections) for(auto const& ip:element)
        { result.push_back(ip.plastic_deformation[0]);result.push_back(ip.plastic_deformation[1]); }
        return result;
    }
    void restoreCheckpointState(std::vector<double> const& state) override
    {
        if(!nonlinear_ && state.empty()) { return; }
        auto const n=accepted_.rotations.size();auto const m=accepted_.sections.size();
        if(state.size()!=3+3*n+4*m || state[0]!=1. || state[1]!=n || state[2]!=m ||
            !std::all_of(state.begin(),state.end(),[](double v){return std::isfinite(v);}))
        { OGS_FATAL("Invalid plate checkpoint state."); }
        std::size_t k=3;
        for(Eigen::Index i=0;i<2*n;++i) { accepted_.displacement[i]=state[k++]; }
        for(Eigen::Index i=0;i<n;++i) { accepted_.rotations[i]=state[k++]; }
        for(auto& element:accepted_.sections) for(auto& ip:element)
        { ip.plastic_deformation[0]=state[k++];ip.plastic_deformation[1]=state[k++]; }
    }
    void integrate(double t,GlobalVector const& x,GlobalVector& b,GlobalMatrix* jac) const override
    {
        auto const values=x.get(indices_);
        Eigen::Map<Eigen::VectorXd const> u(values.data(),values.size());
        if(nonlinear_)
        {
            auto const response=nonlinearTrial(u);
            Eigen::VectorXd const rhs=-response.force;b.add(indices_,rhs);
            if(jac) { jac->add({indices_,indices_},response.tangent); }
        }
        else
        {
            Eigen::VectorXd const rhs=-frame_.stiffness*u+loadFactor(t)*condensed_load_;
            b.add(indices_,rhs);
            if (jac) { jac->add({indices_,indices_},frame_.stiffness); }
        }
    }
    void postAcceptedTimestep(double t,GlobalVector const& x) override
    {
        auto const values=x.get(indices_);
        Eigen::Map<Eigen::VectorXd const> u(values.data(),values.size());
        Eigen::VectorXd rotations=frame_.rotation_from_translation*u+loadFactor(t)*rotation_offset_;
        std::vector<Eigen::Vector<double,6>> nonlinear_forces;
        if(nonlinear_)
        {
            auto response=nonlinearTrial(u);rotations=response.state.rotations;
            nonlinear_forces=std::move(response.local_end_forces);
            accepted_=std::move(response.state);
            for(std::size_t e=0;e<accepted_.sections.size();++e)
            {
                for(int ip=0;ip<2;++ip)
                {
                    auto const& ep=accepted_.sections[e][ip].plastic_deformation;
                    INFO("Accepted plate section {} element {} ip {}: t={}, plastic_strain={}, plastic_curvature={}",mesh_.getName(),e,ip,t,ep[0],ep[1]);
                }
            }
        }
        for (auto const* e:mesh_.getElements())
        {
            Eigen::Vector<double,6> q;
            for(int a=0;a<2;++a)
            {
                auto const i=e->getNode(a)->getID();
                q.segment<2>(3*a)=u.segment<2>(2*i);q[3*a+2]=rotations[i];
            }
            auto const delta=(e->getNode(1)->asEigenVector3d()-e->getNode(0)->asEigenVector3d()).head<2>().eval();
            StructuralElements::Beam2DLocalAssembler const local(delta,data_);
            Eigen::Vector<double,6> const f=(nonlinear_?nonlinear_forces.at(e->getID()):local.localForces(q))-loadFactor(t)*member_loads_.at(e->getID());
            INFO("Accepted plate {} element {}: t={}, N1={}, V1={}, M1={}, N2={}, V2={}, M2={}, theta1={}, theta2={}",mesh_.getName(),e->getID(),t,f[0],f[1],f[2],f[3],f[4],f[5],q[2],q[5]);
        }
    }
private:
    double loadFactor(double t) const { return load_ramp_?load_ramp_->factor(t):1.; }
    StructuralElements::FrameNetworkResponse nonlinearTrial(Eigen::VectorXd const& u) const
    {
        try
        {
            Eigen::VectorXd const guess=accepted_.rotations+frame_.rotation_from_translation*(u-accepted_.displacement);
            return nonlinear_->trial(u,accepted_,guess);
        }
        catch(std::runtime_error const& error)
        { throw NumLib::AssemblyException(error.what()); }
    }
    std::unique_ptr<StructuralElements::NonlinearFrameNetwork> nonlinear_;
    StructuralElements::FrameNetworkState accepted_;
    MeshLib::Mesh const& mesh_;
    StructuralElements::Beam2DProcessData data_;
    std::optional<StructuralElements::DistributedLoadRamp> load_ramp_;
    StructuralElements::CondensedFrame frame_;
    Eigen::VectorXd rotation_offset_,condensed_load_;
    std::vector<Eigen::Vector<double,6>> member_loads_;
    std::vector<GlobalIndexType> indices_;
};
}
