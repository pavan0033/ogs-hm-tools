// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <cmath>
#include <map>
#include "BoundaryCondition.h"
#include "BaseLib/Logging.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "NumLib/TimeStepping/Time.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
/// Saturated Darcy boundary: prescribed submerged pressure; otherwise
/// p<=0, outward mass flow>=0 and p*outward_flow=0. No retention law.
class ReservoirSeepageBoundaryCondition final : public BoundaryCondition
{
public:
    ReservoirSeepageBoundaryCondition(ParameterLib::Parameter<double> const* pressure,
        MeshLib::Mesh const& boundary, MeshLib::Mesh const& bulk,
        NumLib::LocalToGlobalIndexMap const& dof, int variable, int component,
        double start, double end, double pressure_tolerance, double flux_tolerance,
        bool seepage_only=false)
        : _pressure(pressure), _bulk(bulk), _start(start), _end(end),
          _ptol(pressure_tolerance), _qtol(flux_tolerance), _seepage_only(seepage_only)
    {
        if ((!pressure && !seepage_only) || variable!=0 || component!=0 || dof.getNumberOfVariables()!=2 ||
            !std::isfinite(start) || !std::isfinite(end) || end<=start ||
            !std::isfinite(_ptol) || !std::isfinite(_qtol) || _ptol<=0 || _qtol<=0)
        { OGS_FATAL("ReservoirSeepage requires scalar HM pressure, positive tolerances and finite start<end."); }
        auto const* map=boundary.getProperties().getPropertyVector<std::size_t>("bulk_node_ids",MeshLib::MeshItemType::Node,1);
        if (!map || map->size()!=boundary.getNumberOfNodes()) { OGS_FATAL("ReservoirSeepage needs a mapped boundary mesh."); }
        for (auto const* node:boundary.getNodes())
        {
            auto const id=(*map)[node->getID()];
            if (id>=bulk.getNumberOfNodes()) { OGS_FATAL("Invalid seepage boundary map."); }
            auto const index=dof.getGlobalIndex({bulk.getID(),MeshLib::MeshItemType::Node,id},variable,component);
            if (index==NumLib::MeshComponentMap::nop) { continue; }
            if (index<0) { OGS_FATAL("ReservoirSeepage currently requires serial pressure DOFs."); }
            _nodes.push_back(id);_indices.push_back(index);
        }
        _active.assign(_nodes.size(),0.);
    }
    void preTimestep(double t,std::vector<GlobalVector*> const& x,int process_id) override
    {
        if (!enabled(t)) { return; }
        // Warm-start a new phase's face selection from the accepted pressure.
        // A zero-pressure outlet must not become sealed merely because a new
        // phase owns a fresh boundary object.
        for (std::size_t i=0;i<_nodes.size();++i)
        { _active[i]=x[process_id]->get(_indices[i])>=-_ptol?1.:0.; }
    }
    bool hasActiveSet() const override { return true; }
    bool enabled(double t) const { return NumLib::Time{t}>NumLib::Time{_start} && NumLib::Time{t}<=NumLib::Time{_end}; }
    void getEssentialBCValues(double t,GlobalVector const&,NumLib::IndexValueVector<GlobalIndexType>& values) const override
    {
        values.ids.clear();values.values.clear();if (!enabled(t)) { return; }
        for (std::size_t i=0;i<_nodes.size();++i)
        {
            double const p=prescribed(t,i);
            if ((!_seepage_only && p>=0.) || _active[i]!=0.) { values.ids.push_back(_indices[i]);values.values.push_back(std::max(p,0.)); }
        }
    }
    void setPumpingMassRhs(std::map<std::size_t,double> const& rates)
    { _pumping_rhs=rates; }
    bool updateActiveSet(double t,GlobalVector const& x) override
    {
        if (!enabled(t)) { return false; }
        auto const* mass=_bulk.getProperties().getPropertyVector<double>("MassFlowRate",MeshLib::MeshItemType::Node,1);
        if (!mass) { OGS_FATAL("Seepage active set requires fresh HM mass residuals."); }
        bool changed=false;
        for (std::size_t i=0;i<_nodes.size();++i)
        {
            if (!_seepage_only && prescribed(t,i)>=0.) { continue; }
            // MassFlowRate is positive into the domain. Releasing an active
            // zero-pressure node with positive reaction prevents artificial rain.
            bool const was=_active[i]!=0.;
            auto const sink=_pumping_rhs.find(_nodes[i]);
            double const boundary_reaction=(*mass)[_nodes[i]]-
                (sink==_pumping_rhs.end()?0.:sink->second);
            bool const next=was ? !(boundary_reaction>_qtol) : x.get(_indices[i])>_ptol;
            changed|=was!=next;_active[i]=next?1.:0.;
        }
        return changed;
    }
    std::vector<double> checkpointState() const override { return _active; }
    void validateCheckpointState(std::vector<double> const& data,double) const override
    {
        if (data.size()!=_active.size()) { OGS_FATAL("Seepage checkpoint size mismatch."); }
        for (auto v:data) { if (v!=0. && v!=1.) { OGS_FATAL("Invalid seepage checkpoint flag."); } }
    }
    void restoreCheckpointState(std::vector<double> const& data,double t) override
    { validateCheckpointState(data,t);_active=data; }
private:
    double prescribed(double t,std::size_t i) const
    {
        if (_seepage_only) { return 0.; }
        ParameterLib::SpatialPosition pos;pos.setNodeID(_nodes[i]);pos.setCoordinates(*_bulk.getNode(_nodes[i]));
        double const p=(*_pressure)(t,pos).front();if (!std::isfinite(p)) { OGS_FATAL("Nonfinite reservoir pressure."); }return p;
    }
    ParameterLib::Parameter<double> const* _pressure;
    MeshLib::Mesh const& _bulk;
    double _start,_end,_ptol,_qtol;
    bool _seepage_only;
    std::vector<std::size_t> _nodes;
    std::vector<GlobalIndexType> _indices;
    std::vector<double> _active;
    std::map<std::size_t,double> _pumping_rhs;
};
} // namespace ProcessLib
