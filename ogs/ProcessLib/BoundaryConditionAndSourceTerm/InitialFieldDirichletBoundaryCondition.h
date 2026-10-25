// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <cmath>
#include <memory>
#include <optional>
#include <set>
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "BoundaryCondition.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "NumLib/TimeStepping/Time.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
/// Select actual scalar P1 unknowns from the bulk DOF map, without an auxiliary
/// node mesh. Default evaluation at t=0 holds the initial field; a native
/// synchronous initialisation ramp explicitly opts into current-time evaluation.
class InitialFieldDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    InitialFieldDirichletBoundaryCondition(
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bulk,
        NumLib::LocalToGlobalIndexMap const& dof,
        int const variable, int const component,
        std::optional<double> const end_time,
        bool const evaluate_at_current_time = false)
        : _parameter(parameter), _bulk(bulk), _end_time(end_time),
          _evaluate_at_current_time(evaluate_at_current_time)
    {
        if (component != 0 || dof.getNumberOfVariableComponents(variable) != 1 ||
            (end_time && (!std::isfinite(*end_time) || *end_time <= 0)))
        {
            OGS_FATAL("InitialFieldDirichlet requires scalar pressure and a finite positive end_time if supplied.");
        }
        std::set<std::size_t> primary;
        for (auto const* element : bulk.getElements())
        {
            for (unsigned i = 0; i < element->getNumberOfBaseNodes(); ++i)
            {
                primary.insert(element->getNode(i)->getID());
            }
        }
        for (auto const* node : bulk.getNodes())
        {
            auto const index = dof.getGlobalIndex(
                {bulk.getID(), MeshLib::MeshItemType::Node, node->getID()}, variable, component);
            bool const has_dof = index != NumLib::MeshComponentMap::nop;
            if (primary.contains(node->getID()) != has_dof || (has_dof && index < 0))
            {
                OGS_FATAL("InitialFieldDirichlet requires complete serial scalar P1 pressure DOFs on bulk primary nodes.");
            }
            if (has_dof)
            {
                _nodes.push_back(node->getID());
                _indices.push_back(index);
            }
        }
        if (_indices.empty())
        {
            OGS_FATAL("InitialFieldDirichlet selected no pressure DOFs.");
        }
        INFO("InitialFieldDirichlet selected {} primary pressure DOFs internally; no auxiliary node mesh.", _indices.size());
    }

    std::vector<GlobalIndexType> const& indices() const { return _indices; }

    bool requiresExclusiveEssentialDofs() const override { return true; }

    void getEssentialBCValues(
        double const t, GlobalVector const&,
        NumLib::IndexValueVector<GlobalIndexType>& values) const override
    {
        values.ids.clear();
        values.values.clear();
        if (_end_time && NumLib::Time{t} > NumLib::Time{*_end_time})
        {
            return;
        }
        values.ids = _indices;
        ParameterLib::SpatialPosition pos;
        for (auto const node : _nodes)
        {
            pos.setNodeID(node);
            pos.setCoordinates(*_bulk.getNode(node));
            double const value = _parameter(_evaluate_at_current_time ? t : 0., pos).front();
            if (!std::isfinite(value))
            {
                OGS_FATAL("InitialFieldDirichlet initial parameter produced nonfinite pressure.");
            }
            values.values.push_back(value);
        }
    }
private:
    ParameterLib::Parameter<double> const& _parameter;
    MeshLib::Mesh const& _bulk;
    std::optional<double> const _end_time;
    bool const _evaluate_at_current_time;
    std::vector<std::size_t> _nodes;
    std::vector<GlobalIndexType> _indices;
};
}  // namespace ProcessLib
