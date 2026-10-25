// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <vector>

#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace HydroMechanics
{
template <int DisplacementDim>
struct LocalAssemblerInterface : public ProcessLib::LocalAssemblerInterface,
                                 public NumLib::ExtrapolatableElement
{
    virtual std::vector<std::string> captureIntegrationPointCheckpoint() const
    { OGS_FATAL("HM assembler does not expose complete integration-point checkpoints."); }
    virtual void validateIntegrationPointCheckpoint(std::vector<std::string> const&) const
    { OGS_FATAL("HM assembler does not expose checkpoint validation."); }
    virtual void restoreIntegrationPointCheckpoint(std::vector<std::string> const&)
    { OGS_FATAL("HM assembler does not expose complete integration-point checkpoints."); }

    virtual double getInitialisationYieldResidual(double const factor=1.) const
    { OGS_FATAL("This HM assembler does not expose an initialisation gate."); }

    virtual std::size_t setIPDataInitialConditions(
        std::string_view const name, double const* values,
        int const integration_order) = 0;

    virtual std::vector<double> getMaterialStateVariableInternalState(
        typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
            InternalVariable::WriteAccess const& values,
        int const& components) const = 0;

    virtual void acceptMechanicalDrainage(std::size_t mesh_item_id, double t, GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table) = 0;
    virtual std::vector<double> getMechanicalStressBasis() const = 0;
    virtual std::vector<double> getMechanicalPressureValid() const = 0;
    virtual std::vector<double> getMechanicalPressure() const = 0;
    virtual std::vector<double> getAppliedBiotPressure() const = 0;
    virtual std::vector<double> getMechanicalExcess() const = 0;
    virtual std::vector<double> getMechanicalFluidContent() const = 0;
    virtual std::vector<double> getMechanicalReferenceStrain() const = 0;
    virtual std::vector<double> getMechanicalReferencePressure() const = 0;
    virtual std::vector<double> getMechanicalActualReferencePressure() const = 0;
    virtual std::vector<double> getSigma() const = 0;

    virtual std::vector<double> getEpsilon() const = 0;

    virtual std::vector<double> getStrainRateVariable() const = 0;

    virtual std::vector<double> const& getIntPtSigma(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilon(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    // TODO move to NumLib::ExtrapolatableElement
    virtual unsigned getNumberOfIntegrationPoints() const = 0;

    virtual int getMaterialID() const = 0;

    virtual typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned /*integration_point*/) const = 0;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib
