// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <limits>
#include "BoundaryCondition.h"
#include "BaseLib/Logging.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "NumLib/TimeStepping/Time.h"
#include "ProcessLib/HydroMechanics/PumpingWellControl.h"

namespace ProcessLib
{
/// 2D plane-strain, mapped P1 line screen. extraction_rate is total volumetric
/// extraction per out-of-plane metre (m^2/s), not extraction per screen length.
/// The screen rate is lumped by line length. Each node has a local minimum-head
/// constraint; this is not a common-head wellbore/redistribution model.
class PumpingWellBoundaryCondition final : public BoundaryCondition
{
public:
    using State = HydroMechanics::PumpingState;
    PumpingWellBoundaryCondition(MeshLib::Mesh const& screen,
        MeshLib::Mesh const& bulk, NumLib::LocalToGlobalIndexMap const& dof,
        int variable, int component, double extraction, double minimum_head,
        double density, double gravity, double start, double end,
        double pressure_tolerance, double mass_tolerance)
        : _name(screen.getName()), _bulk(bulk), _extraction(extraction),
          _head(minimum_head), _rho(density), _g(gravity), _start(start),
          _end(end), _ptol(pressure_tolerance), _qtol(mass_tolerance)
    {
        if (variable != 0 || component != 0 || dof.getNumberOfVariables() != 2 ||
            bulk.getDimension() != 2 || !std::isfinite(extraction) || extraction<=0 ||
            !std::isfinite(minimum_head) || !std::isfinite(density) || density<=0 ||
            !std::isfinite(gravity) || gravity<=0 || !std::isfinite(start) ||
            !std::isfinite(end) || end<=start || !std::isfinite(_ptol) || _ptol<=0 ||
            !std::isfinite(_qtol) || _qtol<=0)
        { OGS_FATAL("PumpingWell requires serial 2D HM pressure, positive rate/density/gravity/tolerances and finite start<end."); }
        auto const* map=screen.getProperties().getPropertyVector<std::size_t>(
            "bulk_node_ids",MeshLib::MeshItemType::Node,1);
        if (!map || map->size()!=screen.getNumberOfNodes())
        { OGS_FATAL("PumpingWell needs a mapped P1 line screen."); }
        std::map<std::size_t,double> weights;
        double length=0.;
        for (auto const* e:screen.getElements())
        {
            if (e->getDimension()!=1 || e->getNumberOfNodes()!=2)
            { OGS_FATAL("PumpingWell screen must contain only two-node lines."); }
            auto const& a=*e->getNode(0);
            auto const& b=*e->getNode(1);
            double const l=std::hypot(b[0]-a[0],b[1]-a[1],b[2]-a[2]);
            if (!std::isfinite(l) || l<=0)
            { OGS_FATAL("PumpingWell has an invalid screen segment."); }
            length+=l;
            for (unsigned j=0;j<2;++j)
            { weights[(*map)[e->getNode(j)->getID()]]+=l/2.; }
        }
        if (!(length>0.)) { OGS_FATAL("PumpingWell screen is empty."); }
        for (auto const& [id,w]:weights)
        {
            if (id>=bulk.getNumberOfNodes()) { OGS_FATAL("Invalid pumping screen map."); }
            auto const index=dof.getGlobalIndex(
                {bulk.getID(),MeshLib::MeshItemType::Node,id},variable,component);
            if (index==NumLib::MeshComponentMap::nop || index<0)
            { OGS_FATAL("PumpingWell needs serial primary pressure nodes."); }
            _nodes.push_back(id); _indices.push_back(index);
            _rates.push_back(_rho*_extraction*w/length);
            _minimum.push_back(_rho*_g*(_head-(*bulk.getNode(id))[1]));
        }
        _states.assign(_nodes.size(),State::Off);
    }
    std::map<std::size_t,double> massRhs(double t) const
    {
        std::map<std::size_t,double> rhs;
        if (!enabled(t)) { return rhs; }
        for (std::size_t i=0;i<_nodes.size();++i)
        { if (_states[i]==State::Rate) { rhs[_nodes[i]]-=_rates[i]; } }
        return rhs;
    }
    bool overlaps(PumpingWellBoundaryCondition const& other) const
    {
        if (_end<=other._start || other._end<=_start) { return false; }
        for (auto const index:_indices)
        { if (std::find(other._indices.begin(),other._indices.end(),index)!=other._indices.end()) { return true; } }
        return false;
    }
    bool enabled(double t) const
    { return NumLib::Time{t}>NumLib::Time{_start} && NumLib::Time{t}<=NumLib::Time{_end}; }
    bool hasActiveSet() const override { return true; }
    bool requiresExclusiveEssentialDofs() const override { return true; }
    void preTimestep(double t,std::vector<GlobalVector*> const& x,int pid) override
    {
        if (!enabled(t)) { return; }
        for (std::size_t i=0;i<_nodes.size();++i)
        {
            if (_states[i]==State::Off && x[pid]->get(_indices[i])>_minimum[i]+_ptol)
            { _states[i]=State::Rate; }
        }
    }
    void applyNaturalBC(double t,std::vector<GlobalVector*> const&,int,
                        GlobalMatrix*,GlobalVector& b,GlobalMatrix*) override
    {
        if (!enabled(t)) { return; }
        Eigen::VectorXd rhs=Eigen::VectorXd::Zero(_indices.size());
        for (std::size_t i=0;i<_nodes.size();++i)
        { if (_states[i]==State::Rate) { rhs[i]=-_rates[i]; } }
        b.add(_indices,rhs);
    }
    void getEssentialBCValues(double t,GlobalVector const&,
        NumLib::IndexValueVector<GlobalIndexType>& values) const override
    {
        values.ids.clear();values.values.clear();
        if (!enabled(t)) { return; }
        for (std::size_t i=0;i<_nodes.size();++i)
        {
            if (_states[i]==State::MinimumHead)
            { values.ids.push_back(_indices[i]);values.values.push_back(_minimum[i]); }
        }
    }
    bool updateActiveSet(double t,GlobalVector const& x) override
    {
        if (!enabled(t)) { return false; }
        auto const* mass=_bulk.getProperties().getPropertyVector<double>(
            "MassFlowRate",MeshLib::MeshItemType::Node,1);
        if (!mass) { OGS_FATAL("PumpingWell needs fresh HM mass residuals."); }
        bool changed=false;
        for (std::size_t i=0;i<_nodes.size();++i)
        {
            auto const next=HydroMechanics::pumpingState(_states[i],x.get(_indices[i]),
                _minimum[i],(*mass)[_nodes[i]],_rates[i],_ptol,_qtol);
            changed|=next!=_states[i];_states[i]=next;
        }
        return changed;
    }
    void postAcceptedTimestep(double t,std::vector<GlobalVector*> const& x,int pid) override
    {
        if (!enabled(t)) { return; }
        auto const* mass=_bulk.getProperties().getPropertyVector<double>(
            "MassFlowRate",MeshLib::MeshItemType::Node,1);
        if (!mass) { OGS_FATAL("PumpingWell accepted audit needs mass residuals."); }
        double q=0., hmin=std::numeric_limits<double>::infinity(), hmax=-hmin;
        std::size_t limited=0,off=0;
        for (std::size_t i=0;i<_nodes.size();++i)
        {
            if (_states[i]==State::Rate) { q+=_rates[i]/_rho; }
            if (_states[i]==State::MinimumHead) { q-=(*mass)[_nodes[i]]/_rho; }
            double const h=x[pid]->get(_indices[i])/(_rho*_g)+(*_bulk.getNode(_nodes[i]))[1];
            hmin=std::min(hmin,h);hmax=std::max(hmax,h);
            limited+=_states[i]==State::MinimumHead;off+=_states[i]==State::Off;
        }
        INFO("Accepted pumping well {}: t={}, extraction={}, maximum={}, head_min={}, head_max={}, limited={}, off={}",
            _name,t,q,_extraction,hmin,hmax,limited,off);
    }
    std::vector<double> checkpointState() const override
    {
        std::vector<double> v;for(auto s:_states){v.push_back(static_cast<int>(s));}return v;
    }
    void validateCheckpointState(std::vector<double> const& data,double) const override
    {
        if(data.size()!=_states.size()){OGS_FATAL("PumpingWell checkpoint size mismatch.");}
        for(double s:data){if(s!=0.&&s!=1.&&s!=2.){OGS_FATAL("Invalid pumping state.");}}
    }
    void restoreCheckpointState(std::vector<double> const& data,double t) override
    {
        validateCheckpointState(data,t);
        for(std::size_t i=0;i<data.size();++i){_states[i]=static_cast<State>(static_cast<int>(data[i]));}
    }
private:
    std::string _name;
    MeshLib::Mesh const& _bulk;
    double _extraction,_head,_rho,_g,_start,_end,_ptol,_qtol;
    std::vector<std::size_t> _nodes;
    std::vector<GlobalIndexType> _indices;
    std::vector<double> _rates,_minimum;
    std::vector<State> _states;
};
} // namespace ProcessLib
