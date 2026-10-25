// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateAnchorTerm.h"
#include <cmath>

#include "AnchorTerm.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
template <int GlobalDim>
std::unique_ptr<SourceTerm> createAnchorTerm(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
    std::size_t const source_term_mesh_id, const int variable_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Constructing AnchorTerm from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    config.checkConfigParameter("type", "Anchor");

    auto const param_name =
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__Anchor__anchor_force_constant}
        config.getConfigParameter<std::string>("anchor_force_constant");
    DBUG("Using parameter {:s} as anchor stiffness constant.", param_name);

    auto& param = ParameterLib::findParameter<double>(param_name, parameters, 1,
                                                      &st_mesh);

    auto const law = config.getConfigParameterOptional<std::string>("axial_law");
    bool const plastic = law.has_value();
    bool const tension_only = law && *law == "TensionOnly";
    double tension_limit = 0.;
    double compression_limit = 0.;
    if (plastic)
    {
        if (*law != "ElasticPerfectlyPlastic" && !tension_only)
        { OGS_FATAL("Unknown anchor axial_law {}", *law); }
        tension_limit = config.getConfigParameter<double>("tension_limit");
        if (!tension_only)
        { compression_limit = config.getConfigParameter<double>("compression_limit"); }
        if (!std::isfinite(tension_limit) || tension_limit < 0. ||
            !std::isfinite(compression_limit) || compression_limit < 0.)
        { OGS_FATAL("Anchor force limits must be finite and nonnegative."); }
    }

    double const initial_force=config.getConfigParameter<double>("initial_force",0.);
    double const prestress_end_time=config.getConfigParameter<double>("prestress_end_time",-1.);
    double const activation_time=config.getConfigParameter<double>("activation_time",0.);
    if(!std::isfinite(prestress_end_time) || (prestress_end_time!=-1. && prestress_end_time<=activation_time))
    { OGS_FATAL("Anchor prestress_end_time must be finite and greater than activation_time."); }
    if(!std::isfinite(initial_force) || !std::isfinite(activation_time) || activation_time<0. ||
        (plastic && (initial_force>tension_limit || initial_force< -compression_limit)))
    { OGS_FATAL("Anchor prestress must be finite and admissible; activation_time must be finite and nonnegative."); }
    if ((plastic || initial_force!=0. || activation_time>0. || prestress_end_time>=0.) && param.isTimeDependent())
    { OGS_FATAL("History-dependent anchors require time-independent EA."); }

    for (MeshLib::Element const* const element : st_mesh.getElements())
    {
        if (element->getNumberOfNodes() != 2)
        {
            OGS_FATAL(
                "Every anchor element needs to have precisely two nodes.");
        }
        if (plastic && (element->getNode(1)->asEigenVector3d() -
                        element->getNode(0)->asEigenVector3d()).norm() <= 0.)
        { OGS_FATAL("Anchor length must be positive."); }

    }

    return std::make_unique<AnchorTerm<GlobalDim>>(
        std::move(dof_table), source_term_mesh_id, st_mesh, variable_id, param,
        plastic, tension_limit, compression_limit, tension_only,initial_force,activation_time,prestress_end_time);
}

template std::unique_ptr<SourceTerm> createAnchorTerm<2>(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
    std::size_t const source_term_mesh_id, const int variable_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

template std::unique_ptr<SourceTerm> createAnchorTerm<3>(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
    std::size_t const source_term_mesh_id, const int variable_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}  // namespace ProcessLib
