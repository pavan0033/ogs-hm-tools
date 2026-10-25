// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include "Beam2DProcess.h"
#include "BaseLib/ConfigTree.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include <limits>
namespace ProcessLib::StructuralElements
{
inline std::unique_ptr<Process> createBeam2DProcess(std::string const& name,
    MeshLib::Mesh& mesh,std::unique_ptr<AbstractJacobianAssembler>&& jacobian,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type","BEAM_MECHANICS");
    auto const pv_config=config.getConfigSubtree("process_variables");
    auto pv=findProcessVariables(variables,pv_config,{"displacement","rotation"});
    if(pv[0].get().getNumberOfGlobalComponents()!=2 || pv[1].get().getNumberOfGlobalComponents()!=1 ||
       pv[0].get().getShapeFunctionOrder()!=1 || pv[1].get().getShapeFunctionOrder()!=1)
    { OGS_FATAL("BEAM_MECHANICS requires first-order displacement[2] and rotation[1]."); }
    Beam2DProcessData data{config.getConfigParameter<double>("axial_stiffness"),
        config.getConfigParameter<double>("bending_stiffness"),
        config.getConfigParameter<double>("shear_stiffness",std::numeric_limits<double>::infinity())};
    auto const law=config.getConfigParameterOptional<std::string>("section_law");
    if(law)
    {
        if(*law!="DiamondPerfectPlastic")
        { OGS_FATAL("Unknown beam section_law {}",*law); }
        data.axial_capacity=config.getConfigParameter<double>("axial_capacity");
        data.moment_capacity=config.getConfigParameter<double>("moment_capacity");
        if(!std::isfinite(data.axial_capacity) || data.axial_capacity<=0. ||
           !std::isfinite(data.moment_capacity) || data.moment_capacity<=0. ||
           std::isfinite(data.GA))
        { OGS_FATAL("Plastic Euler beams need positive finite capacities and omitted shear_stiffness."); }
    }
    // Reuse the kernel's finite positive-property validation.
    frameStiffness(1.,data.EA,data.EI,data.GA);
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>> pvs;
    pvs.push_back(std::move(pv));
    return std::make_unique<Beam2DProcess>(name,mesh,std::move(jacobian),parameters,
        integration_order,std::move(pvs),data);
}
}
