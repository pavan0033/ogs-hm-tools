// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <algorithm>
#include <array>
#include <Eigen/LU>

#include "SourceTerm.h"
#include "BaseLib/Logging.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/StructuralElements/Contact2D.h"

namespace ProcessLib
{
/// Weighted pair contact. kn/kt are Pa/m, strength Pa, measure m per unit strip.
/// Optional water traction uses the primary P1 soil pressure, not strict-ESA IP
/// excess pressure. It adds no hydraulic conductivity across the interface.
class Interface2D final : public SourceTerm
{
public:
    Interface2D(std::unique_ptr<NumLib::LocalToGlobalIndexMap> map,
                MeshLib::Mesh const& mesh, Eigen::Vector2d normal,
                double kn, double kt, double cohesion, double friction,
                std::string reduction_policy, MeshLib::Mesh const& bulk,
                NumLib::LocalToGlobalIndexMap const& bulk_dof,
                int pressure_variable, unsigned soil_side,
                ParameterLib::Parameter<double> const* initial_normal,
                bool retain_cohesion_when_open=false)
        : SourceTerm(std::move(map)), mesh_(mesh), kn_(kn), kt_(kt),
          cohesion_(cohesion), mu_(friction),
          retain_cohesion_when_open_(retain_cohesion_when_open && cohesion>0.),
          reduction_policy_(std::move(reduction_policy)),
          states_(mesh.getNumberOfElements()),
          preload_(mesh.getNumberOfElements(), 0.),
          active_(mesh.getNumberOfElements(), true)
    {
#ifdef USE_PETSC
        OGS_FATAL("Interface2D currently requires serial Eigen.");
#endif
        if (!normal.allFinite() || std::abs(normal.norm()-1.)>1e-10 || soil_side>1)
        { OGS_FATAL("Interface2D needs a unit normal and soil_side 0 or 1."); }
        Eigen::Vector2d const tangent(-normal[1],normal[0]);
        C_ << normal[0],-normal[0],normal[1],-normal[1],
             -tangent[0],tangent[0],-tangent[1],tangent[1];
        measures_=mesh.getProperties().getPropertyVector<double>("interface_measure");
        if (!measures_ || measures_->size()!=states_.size())
        { OGS_FATAL("Interface2D requires one positive interface_measure per pair."); }
        StructuralElements::contactTrial(Eigen::Vector2d::Zero(),{},kn_,kt_,cohesion_,mu_);
        if (mesh.getProperties().existsPropertyVector<std::size_t>("soil_element_ids"))
        {
            owners_=mesh.getProperties().getPropertyVector<std::size_t>("soil_element_ids");
            if (owners_->size()!=states_.size())
            { OGS_FATAL("Interface soil_element_ids must have one owner per pair."); }
        }
        if (pressure_variable>=0 && !owners_)
        { OGS_FATAL("Wet Interface2D requires soil_element_ids on each pair."); }
        pressure_.resize(pressure_variable>=0 ? states_.size() : 0);
        for (auto const* e : mesh.getElements())
        {
            auto const i=e->getID();
            if (e->getNumberOfNodes()!=2 || !std::isfinite((*measures_)[i]) || (*measures_)[i]<=0.)
            { OGS_FATAL("Interface2D requires two-node pairs and positive finite measures."); }
            if (initial_normal)
            {
                ParameterLib::SpatialPosition pos;
                pos.setCoordinates(*e->getNode(soil_side));
                preload_[i]=(*initial_normal)(0.,pos)[0];
                if (!std::isfinite(preload_[i]) || preload_[i]<0.)
                { OGS_FATAL("Interface initial normal pressure must be finite and non-negative."); }
            }
            if (owners_ && (*owners_)[i]>=bulk.getNumberOfElements())
            { OGS_FATAL("Invalid interface soil owner."); }
            if (pressure_variable<0) { continue; }
            auto const& owner=*bulk.getElement((*owners_)[i]);
            if (owner.getNumberOfBaseNodes()!=3 ||
                (owner.getNumberOfNodes()!=3 && owner.getNumberOfNodes()!=6))
            { OGS_FATAL("Wet Interface2D currently requires T3/T6 soil owners."); }
            auto const a=owner.getNode(0)->asEigenVector3d().head<2>().eval();
            Eigen::Matrix2d J;
            J.col(0)=owner.getNode(1)->asEigenVector3d().head<2>()-a;
            J.col(1)=owner.getNode(2)->asEigenVector3d().head<2>()-a;
            if (std::abs(J.determinant())<1e-18)
            { OGS_FATAL("Degenerate interface soil owner."); }
            if (owner.getNumberOfNodes()==6)
            {
                constexpr std::array<std::array<int,2>,3> edge{{{{0,1}},{{1,2}},{{2,0}}}};
                for (unsigned j=0;j<3;++j)
                {
                    auto const midpoint=((owner.getNode(edge[j][0])->asEigenVector3d()+
                                          owner.getNode(edge[j][1])->asEigenVector3d())/2.).eval();
                    if ((owner.getNode(3+j)->asEigenVector3d()-midpoint).norm()>1e-9*std::max(1.,J.norm()))
                    { OGS_FATAL("Wet Interface2D requires straight-sided T6 geometry."); }
                }
            }
            Eigen::Vector2d const rs=J.fullPivLu().solve(e->getNode(soil_side)->asEigenVector3d().head<2>()-a);
            auto& p=pressure_[i];p.N << 1.-rs.sum(),rs[0],rs[1];
            if (p.N.minCoeff()<-1e-9 || p.N.maxCoeff()>1.+1e-9)
            { OGS_FATAL("Interface soil-side point is outside its declared owner."); }
            for (unsigned j=0;j<3;++j)
            {
                auto const id=bulk_dof.getGlobalIndex({bulk.getID(),MeshLib::MeshItemType::Node,owner.getNode(j)->getID()},pressure_variable,0);
                if (id==NumLib::MeshComponentMap::nop || id<0)
                { OGS_FATAL("Wet interface requires valid serial P1 pressure DOFs."); }
                p.indices.push_back(id);
            }
        }
    }
    bool hasHistory() const override { return true; }
    bool supportsSafety() const override { return !retain_cohesion_when_open_ && reduction_policy_!="Unspecified"; }
    void setStrengthReductionFactor(double factor) override
    { factor_=reduction_policy_=="ReduceWithSoil"?factor:1.; }
    void updateActiveElements(std::vector<std::size_t> const& ids,bool masked) override
    {
        if (!owners_) { return; }
        for (std::size_t i=0;i<active_.size();++i)
        { active_[i]=!masked || std::find(ids.begin(),ids.end(),(*owners_)[i])!=ids.end(); }
    }
    bool supportsCheckpoint() const override { return true; }
    std::vector<double> checkpointState() const override
    {
        std::vector<double> result{retain_cohesion_when_open_?2.:1.};
        if (retain_cohesion_when_open_) { result.push_back(1.); }
        for (auto const& s:states_) { result.push_back(s.plastic_slip); }
        return result;
    }
    void restoreCheckpointState(std::vector<double> const& values) override
    {
        std::size_t const offset=retain_cohesion_when_open_?2:1;
        if (values.size()!=states_.size()+offset ||
            values[0]!=(retain_cohesion_when_open_?2.:1.) ||
            (retain_cohesion_when_open_ && values[1]!=1.) ||
            !std::all_of(values.begin(),values.end(),[](double v){return std::isfinite(v);}))
        { OGS_FATAL("Invalid Interface2D checkpoint state."); }
        for (std::size_t i=0;i<states_.size();++i) { states_[i].plastic_slip=values[i+offset]; }
    }
    void integrate(double,GlobalVector const& x,GlobalVector& b,GlobalMatrix* jac) const override
    {
        for (auto const* e:mesh_.getElements())
        {
            auto const i=e->getID();if (!active_[i]) { continue; }
            auto const ids=NumLib::getIndices(i,*_source_term_dof_table);
            auto const v=x.get(ids);Eigen::Map<Eigen::Vector4d const> u(v.data());
            Eigen::Vector2d deformation=C_*u;deformation[0]+=preload_[i]/kn_;
            auto response=StructuralElements::contactTrial(deformation,states_[i],kn_,kt_,cohesion_/factor_,mu_/factor_,retain_cohesion_when_open_);
            auto const p=waterPressure(i,x);response.traction[0]+=std::max(p,0.);
            double const w=(*measures_)[i];
            Eigen::Vector4d const rhs=-w*C_.transpose()*response.traction;b.add(ids,rhs);
            if (jac)
            {
                Eigen::Matrix4d const J=w*C_.transpose()*response.tangent*C_;jac->add({ids,ids},J);
                if (!pressure_.empty())
                {
                    Eigen::Matrix<double,4,3> const up=w*(p>0.?1.:0.)*C_.row(0).transpose()*pressure_[i].N.transpose();
                    jac->add({ids,pressure_[i].indices},up);
                }
            }
        }
    }
    void postAcceptedTimestep(double t,GlobalVector const& x) override
    {
        for (auto const* e:mesh_.getElements())
        {
            auto const i=e->getID();if (!active_[i]) { continue; }
            auto const ids=NumLib::getIndices(i,*_source_term_dof_table);
            auto const v=x.get(ids);Eigen::Map<Eigen::Vector4d const> u(v.data());
            Eigen::Vector2d const d=C_*u;auto effective=d;effective[0]+=preload_[i]/kn_;
            auto const response=StructuralElements::contactTrial(effective,states_[i],kn_,kt_,cohesion_/factor_,mu_/factor_,retain_cohesion_when_open_);
            states_[i]=response.state;
            INFO("Accepted interface {} pair {}: t={}, closure={}, slip={}, normal={}, shear={}, plastic_slip={}",mesh_.getName(),i,t,d[0],d[1],response.traction[0],response.traction[1],response.state.plastic_slip);
            if (!pressure_.empty())
            { INFO("Accepted interface water {} pair {}: t={}, pressure={}",mesh_.getName(),i,t,std::max(waterPressure(i,x),0.)); }
        }
    }
private:
    struct PressureTrace { Eigen::Vector3d N;std::vector<GlobalIndexType> indices; };
    double waterPressure(std::size_t i,GlobalVector const& x) const
    {
        if (pressure_.empty()) { return 0.; }
        auto const values=x.get(pressure_[i].indices);
        return pressure_[i].N.dot(Eigen::Map<Eigen::Vector3d const>(values.data()));
    }
    MeshLib::Mesh const& mesh_;
    double kn_,kt_,cohesion_,mu_;
    bool retain_cohesion_when_open_;
    std::string reduction_policy_;
    double factor_=1.;
    Eigen::Matrix<double,2,4> C_;
    MeshLib::PropertyVector<double> const* measures_;
    MeshLib::PropertyVector<std::size_t> const* owners_=nullptr;
    std::vector<StructuralElements::ContactState> states_;
    std::vector<double> preload_;
    std::vector<bool> active_;
    std::vector<PressureTrace> pressure_;
};
}
