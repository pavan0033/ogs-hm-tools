// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include "SourceTerm.h"
#include "BaseLib/Logging.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Exceptions.h"
#include "ProcessLib/StructuralElements/EmbeddedPileNetwork.h"
#include <algorithm>
#include <map>
namespace ProcessLib
{
/// Straight independent beam embedded through validated T3/T6 interpolation.
/// Individual pile forces/stiffnesses are divided by spacing only on transfer
/// to the unit-width soil continuum, avoiding double spacing conversion.
class EmbeddedPile2D final:public SourceTermBase
{
public:
    EmbeddedPile2D(MeshLib::Mesh const& soil,MeshLib::Mesh const& pile,
                    NumLib::LocalToGlobalIndexMap const& dofs,int variable,
                    StructuralElements::PileProperties props,std::size_t head,std::size_t tip,
                    double spacing,double activation,std::string reduction_policy,
                    ParameterLib::Parameter<double> const* force_a,
                    ParameterLib::Parameter<double> const* force_l,
                    ParameterLib::Parameter<double> const* displacement_a,
                    ParameterLib::Parameter<double> const* displacement_l,
                    bool bond_head=false,bool use_interaction=false,bool use_stiffness=false,
                    double prestress_force=0.,double prestress_end=-1.,bool prestress_ramp=false,
                    std::size_t external_head=std::numeric_limits<std::size_t>::max())
        :soil_(soil),dofs_(dofs),variable_(variable),has_base_(props.base_stiffness>0. && props.base_limit>0.),prestress_ramp_(prestress_ramp),prestress_force_(prestress_force),prestress_end_(prestress_end),bond_head_(bond_head),mesh_(pile),head_(head),tip_(tip),spacing_(spacing),activation_(activation),reduction_policy_(std::move(reduction_policy)),
         force_a_(force_a),force_l_(force_l),displacement_a_(displacement_a),displacement_l_(displacement_l),external_head_(external_head)
    {
#ifdef USE_PETSC
        OGS_FATAL("EmbeddedPile2D currently requires serial Eigen.");
#endif
        if(!std::isfinite(prestress_force_) || prestress_force_<0. || !std::isfinite(prestress_end_) ||
           (prestress_force_!=0. && prestress_end_<0.) ||
           (prestress_end_!=-1. && (prestress_end_<=activation || !bond_head)))
        { OGS_FATAL("Embedded prestress needs nonnegative finite force, BondedToSoil head and end time after activation."); }
        auto const n=pile.getNumberOfNodes();
        if(external_head_!=std::numeric_limits<std::size_t>::max() &&
           (!bond_head || external_head_>=soil.getNumberOfNodes() || head>=n ||
            (soil.getNode(external_head_)->asEigenVector3d()-pile.getNode(head)->asEigenVector3d()).norm()>1e-8))
        { OGS_FATAL("Embedded explicit head node must coincide with the bonded member head."); }
        if(n<2 || n>200 || head>=n || tip>=n || head==tip || !std::isfinite(spacing) || spacing<=0. ||
           !std::isfinite(activation) || activation<0. || props.lateral_stiffness<=0. ||
           (force_a && displacement_a) || (force_l && displacement_l) ||
           (bond_head && (displacement_a || displacement_l)))
        { OGS_FATAL("Invalid embedded pile size, spacing, activation, lateral stiffness or head controls."); }
        Eigen::Vector2d const delta=(pile.getNode(tip)->asEigenVector3d()-pile.getNode(head)->asEigenVector3d()).head<2>();
        double const length=delta.norm();
        if(length<=0.) { OGS_FATAL("Embedded pile head/tip must be distinct."); }
        axis_=delta/length;lateral_=Eigen::Vector2d(-axis_[1],axis_[0]);
        std::vector<double> s(n);std::vector<int> degree(n,0);
        for(std::size_t i=0;i<n;++i)
        {
            Eigen::Vector2d const d=(pile.getNode(i)->asEigenVector3d()-pile.getNode(head)->asEigenVector3d()).head<2>();
            s[i]=axis_.dot(d);
            if(std::abs(lateral_.dot(d))>1e-8*std::max(1.,length) || s[i]<-1e-10 || s[i]>length+1e-10)
            { OGS_FATAL("EmbeddedPile2D currently requires a straight head-to-tip chain."); }
        }
        std::vector<std::array<std::size_t,2>> members;double total_length=0.;
        for(auto const* e:pile.getElements())
        {
            if(e->getNumberOfNodes()!=2) { OGS_FATAL("Embedded pile members require two nodes."); }
            auto const a=e->getNode(0)->getID(),b=e->getNode(1)->getID();
            members.push_back({a,b});++degree[a];++degree[b];total_length+=std::abs(s[b]-s[a]);
        }
        if(degree[head]!=1 || degree[tip]!=1 || members.size()!=n-1 || std::abs(total_length-length)>1e-8*std::max(1.,length))
        { OGS_FATAL("Embedded pile connectivity must form one straight connected chain."); }
        for(std::size_t i=0;i<n;++i)if(i!=head && i!=tip && degree[i]!=2) { OGS_FATAL("Embedded pile has a branch or unused node."); }
        auto profile=[&](char const* name)
        {
            auto const* values=pile.getProperties().getPropertyVector<double>(name);
            if(!values || values->size()!=members.size() || values->getMeshItemType()!=MeshLib::MeshItemType::Cell || values->getNumberOfGlobalComponents()!=1)
            { OGS_FATAL("Embedded pile profile {} must be scalar cell data.",name); }
            return std::vector<double>(values->begin(),values->end());
        };
        std::vector<double> factors;std::vector<std::array<double,2>> stiffness;
        if(use_interaction) { factors=profile("interaction_factor"); }
        if(use_stiffness)
        {
            auto const EA=profile("pile_axial_stiffness"),EI=profile("pile_bending_stiffness");
            for(std::size_t i=0;i<members.size();++i) { stiffness.push_back({EA[i],EI[i]}); }
        }
        for(std::size_t i=0;i<members.size();++i) { interacting_members_.push_back(factors.empty() || factors[i]>0.); }
        network_=std::make_unique<StructuralElements::EmbeddedPileNetwork>(s,members,head,tip,props,std::move(factors),std::move(stiffness));
        if(external_head_!=std::numeric_limits<std::size_t>::max() && network_->weights()[head_]!=0.)
        { OGS_FATAL("Explicit structural head currently requires a non-interacting free length."); }
        if(prestress_end_>=0. && network_->prestressLength()<=0.)
        { OGS_FATAL("Embedded prestress requires a free length with zero interaction_factor."); }
        accepted_=network_->initialState();
        auto const* element_ids=pile.getProperties().getPropertyVector<std::size_t>("bulk_element_ids");
        auto const* natural=pile.getProperties().getPropertyVector<double>("natural_coordinates");
        if(!element_ids || element_ids->size()!=n || !natural || natural->size()!=3*n)
        { OGS_FATAL("Embedded pile needs nodal bulk_element_ids and three natural_coordinates per node."); }
        struct PointMap { Eigen::VectorXd N;std::array<std::vector<GlobalIndexType>,2> indices; };
        std::vector<PointMap> maps;
        for(std::size_t i=0;i<n;++i)
        {
            if(i==head_ && external_head_!=std::numeric_limits<std::size_t>::max())
            {
                PointMap map;map.N=Eigen::VectorXd::Ones(1);
                for(int c=0;c<2;++c)
                {
                    auto const id=dofs.getGlobalIndex({soil.getID(),MeshLib::MeshItemType::Node,external_head_},variable,c);
                    if(id<0 || id==NumLib::MeshComponentMap::nop) { OGS_FATAL("Embedded explicit head lacks displacement DOFs."); }
                    map.indices[c]={id};indices_.push_back(id);
                }
                maps.push_back(std::move(map));continue;
            }
            if((*element_ids)[i]>=soil.getNumberOfElements()) { OGS_FATAL("Embedded mapping host element out of range."); }
            auto const& host=*soil.getElement((*element_ids)[i]);auto const nodes=host.getNumberOfNodes();
            double const r=(*natural)[3*i],t=(*natural)[3*i+1];
            if(host.getDimension()!=2 || (nodes!=3 && nodes!=6) || !std::isfinite(r) || !std::isfinite(t) || r< -1e-10 || t< -1e-10 || r+t>1+1e-10)
            { OGS_FATAL("Embedded mapping requires valid T3/T6 natural coordinates."); }
            double const a=1-r-t;Eigen::VectorXd N(nodes);
            if(nodes==3) { N<<a,r,t; }
            else { N<<a*(2*a-1),r*(2*r-1),t*(2*t-1),4*a*r,4*r*t,4*t*a; }
            Eigen::Vector3d location=Eigen::Vector3d::Zero();
            for(unsigned j=0;j<nodes;++j) { location+=N[j]*host.getNode(j)->asEigenVector3d(); }
            if((location-pile.getNode(i)->asEigenVector3d()).norm()>1e-8*std::max(1.,length))
            { OGS_FATAL("Embedded natural coordinates do not reproduce the pile point."); }
            PointMap map;map.N=N;
            for(int c=0;c<2;++c)
            {
                map.indices[c]=NumLib::getIndices(variable,c,host.getID(),dofs);
                for(auto const id:map.indices[c])if(id<0 || id==NumLib::MeshComponentMap::nop)
                { OGS_FATAL("Embedded pile requires complete active serial displacement DOFs."); }
                indices_.insert(indices_.end(),map.indices[c].begin(),map.indices[c].end());
            }
            maps.push_back(std::move(map));
        }
        std::sort(indices_.begin(),indices_.end());indices_.erase(std::unique(indices_.begin(),indices_.end()),indices_.end());
        std::map<GlobalIndexType,std::size_t> columns;
        for(std::size_t i=0;i<indices_.size();++i) { columns[indices_[i]]=i; }
        G_=Eigen::MatrixXd::Zero(3*n,indices_.size());
        for(std::size_t i=0;i<n;++i)for(int c=0;c<2;++c)for(Eigen::Index j=0;j<maps[i].N.size();++j)
        {
            auto const column=columns.at(maps[i].indices[c][j]);
            G_(3*i,column)+=axis_[c]*maps[i].N[j];G_(3*i+1,column)+=lateral_[c]*maps[i].N[j];
        }
        original_G_=G_;original_indices_=indices_;
    }
    void updateActiveElements(std::vector<std::size_t> const& ids,bool masked) override
    {
        if(!masked)
        {
            if(mapping_masked_) { G_=original_G_;indices_=original_indices_;active_ids_.clear();mapping_masked_=false; }
            return;
        }
        if(mapping_masked_ && ids==active_ids_) { return; }
        std::vector<bool> active(soil_.getNumberOfElements(),false);
        for(auto id:ids) { if(id>=active.size()) { OGS_FATAL("Embedded active host ID out of range."); }active[id]=true; }
        auto const n=mesh_.getNumberOfNodes();
        auto const* original=mesh_.getProperties().getPropertyVector<std::size_t>("bulk_element_ids");
        auto const* natural=mesh_.getProperties().getPropertyVector<double>("natural_coordinates");
        struct Map { Eigen::VectorXd N;std::array<std::vector<GlobalIndexType>,2> indices; };
        std::vector<Map> maps(n);std::vector<bool> critical(n,false);std::vector<double> coordinate(n);
        std::vector<GlobalIndexType> new_indices;
        for(std::size_t i=0;i<n;++i)
        {
            coordinate[i]=axis_.dot((mesh_.getNode(i)->asEigenVector3d()-mesh_.getNode(head_)->asEigenVector3d()).head<2>());
            critical[i]=network_->weights()[i]>0. || (bond_head_ && i==head_) || (has_base_ && i==tip_);
            if(!critical[i]) { continue; }
            if(i==head_ && external_head_!=std::numeric_limits<std::size_t>::max())
            {
                maps[i].N=Eigen::VectorXd::Ones(1);
                for(int c=0;c<2;++c)
                {
                    auto const id=dofs_.getGlobalIndex({soil_.getID(),MeshLib::MeshItemType::Node,external_head_},variable_,c);
                    if(id<0 || id==NumLib::MeshComponentMap::nop) { OGS_FATAL("Embedded explicit head lacks displacement DOFs."); }
                    maps[i].indices[c]={id};new_indices.push_back(id);
                }
                continue;
            }
            Eigen::Vector2d const point=mesh_.getNode(i)->asEigenVector3d().head<2>();
            std::vector<std::size_t> candidates;
            if((*original)[i]<active.size() && active[(*original)[i]]) { candidates.push_back((*original)[i]); }
            candidates.insert(candidates.end(),ids.begin(),ids.end());bool found=false;
            for(auto eid:candidates)
            {
                auto const& e=*soil_.getElement(eid);auto const count=e.getNumberOfNodes();
                if(e.getDimension()!=2 || (count!=3 && count!=6)) { continue; }
                Eigen::Vector2d rs;
                if(eid==(*original)[i]) { rs<<(*natural)[3*i],(*natural)[3*i+1]; }
                else
                {
                    Eigen::Matrix2d J;J.col(0)=(e.getNode(1)->asEigenVector3d()-e.getNode(0)->asEigenVector3d()).head<2>();
                    J.col(1)=(e.getNode(2)->asEigenVector3d()-e.getNode(0)->asEigenVector3d()).head<2>();
                    if(std::abs(J.determinant())<1e-20) { continue; }
                    rs=J.inverse()*(point-e.getNode(0)->asEigenVector3d().head<2>());
                }
                if(!rs.allFinite() || rs.minCoeff()< -1e-10 || rs.sum()>1.+1e-10) { continue; }
                double const a=1-rs.sum(),r=rs[0],t=rs[1];Eigen::VectorXd N(count);
                if(count==3) { N<<a,r,t; }else { N<<a*(2*a-1),r*(2*r-1),t*(2*t-1),4*a*r,4*r*t,4*t*a; }
                Eigen::Vector2d location=Eigen::Vector2d::Zero();
                for(unsigned j=0;j<count;++j) { location+=N[j]*e.getNode(j)->asEigenVector3d().head<2>(); }
                if((location-point).norm()>1e-8*std::max(1.,(mesh_.getNode(tip_)->asEigenVector3d()-mesh_.getNode(head_)->asEigenVector3d()).norm())) { continue; }
                maps[i].N=N;
                for(int c=0;c<2;++c)
                {
                    maps[i].indices[c]=NumLib::getIndices(variable_,c,eid,dofs_);
                    for(auto index:maps[i].indices[c])
                    {
                        if(index<0 || index==NumLib::MeshComponentMap::nop) { OGS_FATAL("Embedded active mapping lacks displacement DOFs."); }
                        new_indices.push_back(index);
                    }
                }
                if(eid!=(*original)[i]) { INFO("Embedded pile {} node {} remapped from host {} to active host {}",mesh_.getName(),i,(*original)[i],eid); }
                found=true;break;
            }
            if(!found) { OGS_FATAL("Embedded interacting or bonded point {} has no active soil host; excavation removed required support.",i); }
        }
        validateInteractingCoverage(ids);
        std::sort(new_indices.begin(),new_indices.end());new_indices.erase(std::unique(new_indices.begin(),new_indices.end()),new_indices.end());
        std::map<GlobalIndexType,std::size_t> columns;
        for(std::size_t j=0;j<new_indices.size();++j) { columns[new_indices[j]]=j; }
        Eigen::MatrixXd interpolation=Eigen::MatrixXd::Zero(3*n,new_indices.size());
        for(std::size_t i=0;i<n;++i)if(critical[i])
        {
            for(int c=0;c<2;++c)for(Eigen::Index j=0;j<maps[i].N.size();++j)
            {
                auto const column=columns.at(maps[i].indices[c][j]);
                interpolation(3*i,column)+=axis_[c]*maps[i].N[j];interpolation(3*i+1,column)+=lateral_[c]*maps[i].N[j];
            }
        }
        // A free rod has no soil interaction. Use its endpoint chord for any
        // installation reference rather than displacement DOFs inside a void.
        for(std::size_t i=0;i<n;++i)if(!critical[i])
        {
            std::optional<std::size_t> left,right;
            for(std::size_t j=0;j<n;++j)if(critical[j])
            {
                if(coordinate[j]<=coordinate[i] && (!left || coordinate[j]>coordinate[*left])) { left=j; }
                if(coordinate[j]>=coordinate[i] && (!right || coordinate[j]<coordinate[*right])) { right=j; }
            }
            if(!left) { left=right; }if(!right) { right=left; }
            if(!left || !right) { OGS_FATAL("Embedded free length lacks active reference endpoints."); }
            double const fraction=*left==*right?0.:(coordinate[i]-coordinate[*left])/(coordinate[*right]-coordinate[*left]);
            interpolation.block(3*i,0,2,interpolation.cols())=(1-fraction)*interpolation.block(3*(*left),0,2,interpolation.cols())+fraction*interpolation.block(3*(*right),0,2,interpolation.cols());
        }
        G_=std::move(interpolation);indices_=std::move(new_indices);active_ids_=ids;mapping_masked_=true;
    }
    bool hasHistory() const override { return true; }
    bool supportsSafety() const override { return reduction_policy_!="Unspecified"; }
    void setStrengthReductionFactor(double factor) override
    {
        if(factor>1. && prestress_end_>=0. && last_time_<prestress_end_-1e-12*std::max(1.,prestress_end_))
        { OGS_FATAL("Complete embedded prestress lock-off before Safety."); }
        factor_=reduction_policy_=="ReduceWithSoil"?factor:1.;
    }
    bool supportsCheckpoint() const override { return true; }
    std::vector<double> checkpointState() const override
    {
        std::vector<double> out{2.,static_cast<double>(accepted_.shaft.size()),last_time_};
        out.insert(out.end(),accepted_.q.data(),accepted_.q.data()+accepted_.q.size());
        out.insert(out.end(),accepted_.reference.data(),accepted_.reference.data()+accepted_.reference.size());
        for(auto const& s:accepted_.shaft) { out.push_back(s.plastic_slip); }
        for(auto const& s:accepted_.lateral) { out.push_back(s.plastic_slip); }
        out.push_back(accepted_.base.plastic_slip);out.push_back(accepted_.prestress_strain);return out;
    }
    void restoreCheckpointState(std::vector<double> const& values) override
    {
        auto const n=accepted_.shaft.size();
        if((values.size()!=4+8*n && values.size()!=5+8*n) ||
            (values[0]!=1. && values[0]!=2.) ||
            (values[0]==1. && values.size()!=4+8*n) || (values[0]==2. && values.size()!=5+8*n) ||
            values[1]!=n || values[2]<0. ||
            !std::all_of(values.begin(),values.end(),[](double v){return std::isfinite(v);}))
        { OGS_FATAL("Invalid embedded pile checkpoint state."); }
        last_time_=values[2];std::size_t k=3;
        for(Eigen::Index i=0;i<accepted_.q.size();++i) { accepted_.q[i]=values[k++]; }
        for(Eigen::Index i=0;i<accepted_.reference.size();++i) { accepted_.reference[i]=values[k++]; }
        for(auto& s:accepted_.shaft) { s.plastic_slip=values[k++]; }
        for(auto& s:accepted_.lateral) { s.plastic_slip=values[k++]; }
        accepted_.base.plastic_slip=values[k++];
        accepted_.prestress_strain=values[0]==2.?values[k]:0.;
    }
    void integrate(double t,GlobalVector const& x,GlobalVector& b,GlobalMatrix* jac) const override
    {
        if(inactive(t)) { return; }
        if(activation_>0. && last_time_<activation_-1e-12*std::max(1.,activation_))
        { OGS_FATAL("Embedded activation_time must coincide with an accepted boundary."); }
        if(prestress_end_>=0. && t>prestress_end_+1e-12*std::max(1.,prestress_end_) &&
           last_time_<prestress_end_-1e-12*std::max(1.,prestress_end_))
        { OGS_FATAL("Embedded prestress_end_time must coincide with an accepted boundary."); }
        auto const v=x.get(indices_);Eigen::Map<Eigen::VectorXd const> u(v.data(),v.size());
        auto const response=trial(t,G_*u);
        Eigen::VectorXd const rhs=G_.transpose()*response.interface_force/spacing_;b.add(indices_,rhs);
        if(jac) { Eigen::MatrixXd const J=G_.transpose()*response.soil_tangent*G_/spacing_;jac->add({indices_,indices_},J); }
    }
    void postAcceptedTimestep(double t,GlobalVector const& x) override
    {
        auto const v=x.get(indices_);Eigen::Map<Eigen::VectorXd const> u(v.data(),v.size());
        if(inactive(t)) { accepted_.q=G_*u;accepted_.reference=accepted_.q;last_time_=t;return; }
        auto response=trial(t,G_*u);accepted_=std::move(response.state);last_time_=t;
        auto const control=controls(t);
        INFO("Accepted pile {}: t={}, head_axial={}, head_lateral={}, axial_force={}, lateral_force={}, base_force={}, shaft_sum={}",mesh_.getName(),t,accepted_.q[3*head_],accepted_.q[3*head_+1],response.head_reaction[3*head_]+control.external_force[3*head_],response.head_reaction[3*head_+1]+control.external_force[3*head_+1],response.base_force,response.shaft_force.sum());
        for(std::size_t i=0;i<accepted_.shaft.size();++i)
        { INFO("Accepted pile node {} {}: t={}, axial={}, lateral={}, shaft_force={}, lateral_force={}, shaft_plastic_slip={}, lateral_plastic_slip={}",mesh_.getName(),i,t,accepted_.q[3*i],accepted_.q[3*i+1],response.shaft_force[i],response.lateral_force[i],accepted_.shaft[i].plastic_slip,accepted_.lateral[i].plastic_slip); }
    }
private:
    void validateInteractingCoverage(std::vector<std::size_t> const& ids) const
    {
        for(auto const* member:mesh_.getElements())
        {
            if(!interacting_members_.at(member->getID())) { continue; }
            Eigen::Vector2d const start=member->getNode(0)->asEigenVector3d().head<2>();
            Eigen::Vector2d const delta=(member->getNode(1)->asEigenVector3d()-member->getNode(0)->asEigenVector3d()).head<2>();
            std::vector<std::pair<double,double>> intervals;
            for(auto eid:ids)
            {
                auto const& e=*soil_.getElement(eid);auto const count=e.getNumberOfNodes();
                if(e.getDimension()!=2 || (count!=3 && count!=6)) { continue; }
                if(count==6)
                {
                    for(int j=0;j<3;++j)
                    {
                        auto const midpoint=(.5*(e.getNode(j)->asEigenVector3d()+e.getNode((j+1)%3)->asEigenVector3d())).eval();
                        if((midpoint-e.getNode(3+j)->asEigenVector3d()).norm()>1e-8)
                        { OGS_FATAL("Active embedded remapping requires straight-sided T3/T6 host cells."); }
                    }
                }
                Eigen::Matrix2d J;J.col(0)=(e.getNode(1)->asEigenVector3d()-e.getNode(0)->asEigenVector3d()).head<2>();
                J.col(1)=(e.getNode(2)->asEigenVector3d()-e.getNode(0)->asEigenVector3d()).head<2>();
                if(std::abs(J.determinant())<1e-20) { continue; }
                auto const inverse=J.inverse().eval();Eigen::Vector2d const rs=inverse*(start-e.getNode(0)->asEigenVector3d().head<2>()),dr=inverse*delta;
                Eigen::Vector3d bary(1-rs.sum(),rs[0],rs[1]),slope(-dr.sum(),dr[0],dr[1]);
                double low=0.,high=1.;bool valid=true;
                for(int k=0;k<3;++k)
                {
                    if(std::abs(slope[k])<1e-14) { if(bary[k]<-1e-10) { valid=false;break; } }
                    else if(slope[k]>0.) { low=std::max(low,(-1e-10-bary[k])/slope[k]); }
                    else { high=std::min(high,(-1e-10-bary[k])/slope[k]); }
                }
                if(valid && low<=high) { intervals.emplace_back(low,high); }
            }
            std::sort(intervals.begin(),intervals.end());double covered=0.;
            for(auto const& [low,high]:intervals)
            {
                if(low>covered+1e-8) { break; }
                covered=std::max(covered,high);
            }
            if(covered<1.-1e-8)
            { OGS_FATAL("Embedded interacting member {} crosses an excavated region without active soil support.",member->getID()); }
        }
    }
    bool inactive(double t) const { return activation_>0. && t<=activation_+1e-12*std::max(1.,activation_); }
    StructuralElements::PileControls controls(double t) const
    {
        ParameterLib::SpatialPosition pos;pos.setCoordinates(*mesh_.getNode(head_));
        auto value=[&](ParameterLib::Parameter<double> const* p){return p?(*p)(t,pos)[0]:0.;};
        StructuralElements::PileControls out;out.bond_head_to_soil=bond_head_;out.external_force=Eigen::VectorXd::Zero(accepted_.q.size());
        out.external_force[3*head_]=value(force_a_);out.external_force[3*head_+1]=value(force_l_);
        if(prestress_end_>=0. && !inactive(t) && t<=prestress_end_+1e-12*std::max(1.,prestress_end_))
        { out.prestress_force=prestress_force_*(prestress_ramp_?std::clamp((t-activation_)/(prestress_end_-activation_),0.,1.):1.); }
        if(displacement_a_) { out.axial_displacement=value(displacement_a_); }
        if(displacement_l_) { out.lateral_displacement=value(displacement_l_); }
        return out;
    }
    StructuralElements::PileResponse trial(double t,Eigen::VectorXd const& soil) const
    {
        try { return network_->trial(soil,accepted_,controls(t),factor_); }
        catch(std::runtime_error const& error) { throw NumLib::AssemblyException(error.what()); }
    }
    MeshLib::Mesh const& soil_;
    NumLib::LocalToGlobalIndexMap const& dofs_;
    int variable_;
    bool has_base_,mapping_masked_=false;
    std::vector<std::size_t> active_ids_;
    std::vector<bool> interacting_members_;
    std::vector<GlobalIndexType> original_indices_;
    Eigen::MatrixXd original_G_;
    bool prestress_ramp_;
    double prestress_force_,prestress_end_;
    bool bond_head_;
    MeshLib::Mesh const& mesh_;
    std::size_t head_,tip_;
    double spacing_,activation_,last_time_=0.;
    std::string reduction_policy_;
    double factor_=1.;
    Eigen::Vector2d axis_,lateral_;
    Eigen::MatrixXd G_;
    std::vector<GlobalIndexType> indices_;
    std::unique_ptr<StructuralElements::EmbeddedPileNetwork> network_;
    StructuralElements::PileState accepted_;
    ParameterLib::Parameter<double> const* force_a_;
    ParameterLib::Parameter<double> const* force_l_;
    ParameterLib::Parameter<double> const* displacement_a_;
    ParameterLib::Parameter<double> const* displacement_l_;
    std::size_t external_head_;
};
}
