// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Eigenvalues>

#include "HydroMechanicsFEM.h"
#include "HydroMechanicsProcessData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Phase.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/GetFluidDensityAndViscosity.h"
#include "MaterialLib/MPL/Utils/GetSymmetricTensor.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Fem/Interpolation.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"

namespace ProcessLib
{
namespace HydroMechanics
{
namespace MPL = MaterialPropertyLib;

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
HydroMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                             DisplacementDim>::
    HydroMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        HydroMechanicsProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_method),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   _integration_method);

    auto const shape_matrices_p =
        NumLib::initShapeMatrices<ShapeFunctionPressure,
                                  ShapeMatricesTypePressure, DisplacementDim>(
            e, is_axially_symmetric, _integration_method);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            e.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data.integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        // Initialize current time step values
        static const int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
        ip_data.sigma_eff.setZero(kelvin_vector_size);
        ip_data.eps.setZero(kelvin_vector_size);

        // Previous time step values are not initialized and are set later.
        ip_data.eps_prev.resize(kelvin_vector_size);
        ip_data.sigma_eff_prev.resize(kelvin_vector_size);

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_x_prev,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    assert(local_x.size() == pressure_size + displacement_size);
    bool const steady_seepage = _process_data.steady_state_seepage.active(t);
    bool const transient_seepage = _process_data.transient_seepage.active(t);

    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto p_prev =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x_prev.data() + pressure_index,
                                  pressure_size);
    auto u_prev =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x_prev.data() + displacement_index,
                                      displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_Jac_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size>>(
        local_rhs_data, displacement_size + pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType add_p_derivative =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kup = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, pressure_size>::Zero(displacement_size,
                                                    pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        pressure_size, displacement_size>
        Kpu = ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, displacement_size>::Zero(pressure_size,
                                                    displacement_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        pressure_size, displacement_size>
        Kpu_k = ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, displacement_size>::Zero(pressure_size,
                                                    displacement_size);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            _element.getID());

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto const b = (_process_data.initialisationLoadFactor(t) *
                    _process_data.specific_body_force).eval();
    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& solid = medium->phase(MaterialPropertyLib::PhaseName::Solid);
    auto const& fluid = fluidPhase(*medium);
    MPL::VariableArray vars;
    auto& phase_pressure = medium->hasPhase(MaterialPropertyLib::PhaseName::Gas)
                               ? vars.gas_phase_pressure
                               : vars.liquid_phase_pressure;

    auto const T_ref =
        medium->property(MPL::PropertyType::reference_temperature)
            .template value<double>(vars, x_position, t, dt);
    vars.temperature = T_ref;

    auto const& identity2 = Invariants::identity2;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;

        x_position = {
            std::nullopt, _element.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, N_u))};

        auto const x_coord = x_position.getCoordinates().value()[0];

        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;
        eps -= _ip_data[ip].birth_strain;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;

        double const p_int_pt = N_p.dot(p);
        if (_process_data.mechanical_drainage)
        {
            _ip_data[ip].mechanical_pressure = p_int_pt;
            _ip_data[ip].mechanical_excess = 0.;
            _ip_data[ip].mechanical_fluid_content = 0.;
        }

        phase_pressure = p_int_pt;

        auto const C_el = _ip_data[ip].computeElasticTangentStiffness(
            t, x_position, dt, T_ref);
        auto const K_S = solid_material.getBulkModulus(t, x_position, &C_el);

        auto const alpha = usesTotalStress() ? 0. : medium->property(MPL::PropertyType::biot_coefficient)
                               .template value<double>(vars, x_position, t, dt);

        auto const rho_sr =
            solid.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);
        auto const porosity =
            medium->property(MPL::PropertyType::porosity)
                .template value<double>(vars, x_position, t, dt);

        auto [rho_fr, mu] =
            MaterialPropertyLib::getFluidDensityAndViscosity(t, dt, x_position,
                                                             fluid, vars);
        if (_process_data.steady_state_seepage.active(t) ||
            (_process_data.safety && _process_data.safety->active(t) &&
             _process_data.steady_state_seepage.active(_process_data.safety->physical_time)) ||
            (_process_data.mechanical_drainage && _process_data.mechanical_drainage->active(t) &&
             _process_data.steady_state_seepage.active(_process_data.mechanical_drainage->reference_time)))
        {
            // Steady engineering seepage uses one constant reference density.
            // The same EOS still supplies transient density/compressibility.
            auto reference_vars=vars;
            reference_vars.liquid_phase_pressure=0.;
            reference_vars.gas_phase_pressure=0.;
            rho_fr=fluid.property(MPL::PropertyType::density)
                .template value<double>(reference_vars,x_position,t,dt);
        }

        auto const beta_p = steady_seepage ? 0. :
            fluid.property(MPL::PropertyType::density)
                .template dValue<double>(vars, _process_data.phase_variable,
                                         x_position, t, dt) /
            rho_fr;
        double const storage_compressibility = steady_seepage ? 0. : beta_p;

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            auto const sigma_total =
                (_ip_data[ip].sigma_eff - alpha *
                 (_process_data.use_positive_pore_pressure_in_mechanics
                      ? std::max(p_int_pt, 0.) : p_int_pt) * identity2).eval();

            vars.total_stress.emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    sigma_total));
        }
        // For strain dependent permeability
        vars.volumetric_strain = Invariants::trace(eps);
        vars.equivalent_plastic_strain =
            _ip_data[ip].material_state_variables->getEquivalentPlasticStrain();
        vars.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        auto const K = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(vars, x_position, t, dt));
        auto const dkde = MPL::formEigenTensor<
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim)>(
            (*medium)[MPL::PropertyType::permeability].dValue(
                vars, MPL::Variable::mechanical_strain, x_position, t, dt));

        auto const K_over_mu = K / mu;

        if(_process_data.excavation_reference_assembly)
        { _ip_data[ip].sigma_eff=_ip_data[ip].sigma_eff_prev; }
        auto C = _process_data.excavation_reference_assembly ? C_el :
            _ip_data[ip].updateConstitutiveRelation(vars, t, x_position,dt,u,T_ref,
                !_process_data.residual_screening_assembly);

        //
        // displacement equation, displacement part
        //
        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += B.transpose() * C * B * w;

        double const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() -=
            (B.transpose() * sigma_eff - N_u_op(N_u).transpose() * rho * b) * w;

        //
        // displacement equation, pressure part
        //
        if (usesTotalStress())
        {
            // The constitutive stress is total from t=0; no fluid terms.
            continue;
        }
        if (_process_data.use_positive_pore_pressure_in_mechanics &&
            !(_process_data.mechanical_drainage && _process_data.mechanical_drainage->active(t)))
        {
            double const mechanical_pressure = std::max(p_int_pt, 0.);
            _ip_data[ip].mechanical_pressure = mechanical_pressure;
            // Clip the interpolated IP pressure, not the nodal values.
            local_rhs.template segment<displacement_size>(displacement_index)
                .noalias() += B.transpose() * alpha * identity2 * mechanical_pressure * w;
            // Semismooth convention: derivative zero at p=0.
            if (p_int_pt > 0.)
            {
                Kup.noalias() += B.transpose() * alpha * identity2 * N_p * w;
            }
        }
        else if (!_process_data.use_positive_pore_pressure_in_mechanics)
        {
            Kup.noalias() += B.transpose() * alpha * identity2 * N_p * w;
        }

        //
        // pressure equation, pressure part.
        //
        if (_process_data.mechanical_drainage && _process_data.mechanical_drainage->active(t))
        {
            auto const material = _process_data.mechanical_drainage->materialAt(getMaterialID(),t).resolved(alpha, porosity, K_S);
            auto& state = _ip_data[ip].mechanical_drainage;
            auto const response = state.evaluate(p_int_pt, Invariants::trace(eps), alpha, material);
            auto const previous = state.evaluate(N_p.dot(p_prev),
                Invariants::trace((B * u_prev).eval()), alpha, material);
            _ip_data[ip].mechanical_pressure = response.pressure;
            _ip_data[ip].mechanical_excess = response.excess;
            _ip_data[ip].mechanical_fluid_content = response.fluid_content;
            // A/B ignore initial reference suction, not the generated signed
            // increment. This retains full undrained water stiffness on both
            // sides of the phreatic surface. Legacy ClipActual is unchanged.
            bool const positive = _process_data.use_positive_pore_pressure_in_mechanics;
            double const applied = positive ? material.appliedPressure(response.pressure,
                response.excess, state.referencePressure(), true) : response.excess;
            double const derivative = material.pressureDerivative(response.pressure, positive);
            local_rhs.template segment<displacement_size>(displacement_index).noalias() +=
                B.transpose() * alpha * applied * identity2 * w;
            local_Jac.template block<displacement_size, displacement_size>(displacement_index, displacement_index).noalias() +=
                B.transpose() * derivative * response.volumetric_stiffness * identity2 * identity2.transpose() * B * w;
            local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
                rho_fr * N_p.transpose() * (response.fluid_content - previous.fluid_content) / dt * w;
            local_Jac.template block<pressure_size, pressure_size>(pressure_index, pressure_index).noalias() +=
                rho_fr * N_p.transpose() * N_p / (material.biot_modulus * dt) * w;
            if (material.type == MechanicalDrainageType::DrainedEffectiveStress)
            {
                local_Jac.template block<pressure_size, displacement_size>(pressure_index, displacement_index).noalias() +=
                    rho_fr * alpha * N_p.transpose() * identity2.transpose() * B / dt * w;
            }
            continue; // Explicit mechanical idealisation: no Darcy flow.
        }

        laplace_p.noalias() +=
            rho_fr * dNdx_p.transpose() * K_over_mu * dNdx_p * w;

        if (!steady_seepage)
        {
        storage_p.noalias() +=
            rho_fr * N_p.transpose() * N_p * w *
            ((alpha - porosity) * (1.0 - alpha) / K_S + porosity * storage_compressibility);

        }
        // Density dependence of Darcy transport also applies during steady
        // seepage. Storage and volumetric strain rate alone are omitted there.
        // As in stock HM, pressure derivatives of storage are neglected.
        add_p_derivative.noalias() += rho_fr * beta_p * dNdx_p.transpose() *
                                      K_over_mu *
                                      (dNdx_p * p - 2.0 * rho_fr * b) * N_p * w;

        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_fr * rho_fr * K_over_mu * b * w;

        //
        // pressure equation, displacement part.
        //
        // One-way transient seepage retains storage but excludes strain-rate
        // feedback. Darcy transport must not depend on deformation either.
        if (!steady_seepage && !transient_seepage)
        {
        Kpu.noalias() +=
            rho_fr * alpha * N_p.transpose() * identity2.transpose() * B * w;

        Kpu_k.noalias() +=
            dNdx_p.transpose() *
            MathLib::KelvinVector::liftVectorToKelvin<DisplacementDim>(
                dNdx_p * p - rho_fr * b) *
            dkde * B * rho_fr / mu * w;
        }
    }
    // displacement equation, pressure part
    local_Jac
        .template block<displacement_size, pressure_size>(displacement_index,
                                                          pressure_index)
        .noalias() = -Kup;

    if (_process_data.mass_lumping)
    {
        storage_p = storage_p.colwise().sum().eval().asDiagonal();

        if constexpr (pressure_size == displacement_size)
        {
            Kpu = Kpu.colwise().sum().eval().asDiagonal();
            Kpu_k = Kpu_k.colwise().sum().eval().asDiagonal();
        }
    }

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() += laplace_p + storage_p / dt + add_p_derivative;

    // pressure equation, displacement part.
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() += Kpu / dt + Kpu_k;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p + storage_p * (p - p_prev) / dt + Kpu * (u - u_prev) / dt;

    // displacement equation
    if (!_process_data.use_positive_pore_pressure_in_mechanics)
    {
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() += Kup * p;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double> const& HydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, DisplacementDim>::
    getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const
{
    int const hydraulic_process_id = _process_data.hydraulic_process_id;
    auto const indices =
        NumLib::getIndices(_element.getID(), *dof_table[hydraulic_process_id]);
    assert(!indices.empty());
    auto const local_x = x[hydraulic_process_id]->get(indices);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    if (usesTotalStress() || (_process_data.mechanical_drainage && _process_data.mechanical_drainage->active(t)))
    {
        cache.assign(DisplacementDim * n_integration_points, 0.);
        return cache;
    }
    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& fluid = fluidPhase(*medium);
    MPL::VariableArray vars;
    auto& phase_pressure = medium->hasPhase(MaterialPropertyLib::PhaseName::Gas)
                               ? vars.gas_phase_pressure
                               : vars.liquid_phase_pressure;

    // TODO (naumov) Temporary value not used by current material models. Need
    // extension of secondary variables interface.
    double const dt = std::numeric_limits<double>::quiet_NaN();
    vars.temperature =
        medium->property(MPL::PropertyType::reference_temperature)
            .template value<double>(vars, x_position, t, dt);

    auto const& identity2 = Invariants::identity2;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position = {
            std::nullopt, _element.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, _ip_data[ip].N_u))};

        double const p_int_pt = _ip_data[ip].N_p.dot(p);

        phase_pressure = p_int_pt;

        auto const alpha = usesTotalStress() ? 0. : medium->property(MPL::PropertyType::biot_coefficient)
                               .template value<double>(vars, x_position, t, dt);

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        auto const sigma_total =
            (_ip_data[ip].sigma_eff - alpha *
                 (_process_data.use_positive_pore_pressure_in_mechanics
                      ? std::max(p_int_pt, 0.) : p_int_pt) * identity2).eval();
        vars.total_stress.emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_total));

        // For strain dependent permeability
        vars.volumetric_strain = Invariants::trace(_ip_data[ip].eps);
        vars.equivalent_plastic_strain =
            _ip_data[ip].material_state_variables->getEquivalentPlasticStrain();
        vars.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                _ip_data[ip].eps);

        auto const K = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(vars, x_position, t, dt));

        auto [rho_fr, mu] =
            MaterialPropertyLib::getFluidDensityAndViscosity(t, dt, x_position,
                                                             fluid, vars);
        if (_process_data.steady_state_seepage.active(t) ||
            (_process_data.safety && _process_data.safety->active(t) &&
             _process_data.steady_state_seepage.active(_process_data.safety->physical_time)) ||
            (_process_data.mechanical_drainage && _process_data.mechanical_drainage->active(t) &&
             _process_data.steady_state_seepage.active(_process_data.mechanical_drainage->reference_time)))
        {
            // Steady engineering seepage uses one constant reference density.
            // The same EOS still supplies transient density/compressibility.
            auto reference_vars=vars;
            reference_vars.liquid_phase_pressure=0.;
            reference_vars.gas_phase_pressure=0.;
            rho_fr=fluid.property(MPL::PropertyType::density)
                .template value<double>(reference_vars,x_position,t,dt);
        }

        auto const K_over_mu = K / mu;

        auto const b = (_process_data.initialisationLoadFactor(t) *
                    _process_data.specific_body_force).eval();

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -K_over_mu * dNdx_p * p + K_over_mu * rho_fr * b;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForPressureEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data)
{
    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<pressure_size>>(
            local_b_data, pressure_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const p = local_x.template segment<pressure_size>(pressure_index);

    auto const p_prev =
        local_x_prev.template segment<pressure_size>(pressure_index);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, pressure_size>>(local_Jac_data, pressure_size,
                                           pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType add_p_derivative =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            _element.getID());

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& fluid = fluidPhase(*medium);
    MPL::VariableArray vars;
    auto& phase_pressure = medium->hasPhase(MaterialPropertyLib::PhaseName::Gas)
                               ? vars.gas_phase_pressure
                               : vars.liquid_phase_pressure;

    auto const T_ref =
        medium->property(MPL::PropertyType::reference_temperature)
            .template value<double>(vars, x_position, t, dt);
    vars.temperature = T_ref;

    auto const& identity2 = Invariants::identity2;

    auto const staggered_scheme =
        std::get<Staggered>(_process_data.coupling_scheme);
    auto const fixed_stress_stabilization_parameter =
        staggered_scheme.fixed_stress_stabilization_parameter;
    auto const fixed_stress_over_time_step =
        staggered_scheme.fixed_stress_over_time_step;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;

        x_position = {
            std::nullopt, _element.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionPressure,
                                               ShapeMatricesTypePressure>(
                    _element, _ip_data[ip].N_p))};

        double const p_int_pt = N_p.dot(p);

        phase_pressure = p_int_pt;

        auto const C_el = _ip_data[ip].computeElasticTangentStiffness(
            t, x_position, dt, T_ref);
        auto const K_S = solid_material.getBulkModulus(t, x_position, &C_el);

        auto const alpha_b =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, x_position, t, dt);

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        auto const sigma_total =
            (_ip_data[ip].sigma_eff - alpha_b *
                 (_process_data.use_positive_pore_pressure_in_mechanics
                      ? std::max(p_int_pt, 0.) : p_int_pt) * identity2).eval();
        vars.total_stress.emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_total));

        // For strain dependent permeability
        vars.volumetric_strain = Invariants::trace(_ip_data[ip].eps);
        vars.equivalent_plastic_strain =
            _ip_data[ip].material_state_variables->getEquivalentPlasticStrain();
        vars.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                _ip_data[ip].eps);

        auto const K = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(vars, x_position, t, dt));
        auto const porosity =
            medium->property(MPL::PropertyType::porosity)
                .template value<double>(vars, x_position, t, dt);

        auto const [rho_fr, mu] =
            MaterialPropertyLib::getFluidDensityAndViscosity(t, dt, x_position,
                                                             fluid, vars);

        auto const beta_p =
            fluid.property(MPL::PropertyType::density)
                .template dValue<double>(vars, _process_data.phase_variable,
                                         x_position, t, dt) /
            rho_fr;
        double const storage_compressibility = beta_p;

        auto const K_over_mu = K / mu;

        laplace.noalias() +=
            rho_fr * dNdx_p.transpose() * K_over_mu * dNdx_p * w;

        // Artificial compressibility from the fixed stress splitting:
        auto const beta_FS =
            fixed_stress_stabilization_parameter * alpha_b * alpha_b / K_S;

        storage.noalias() += rho_fr * N_p.transpose() * N_p * w *
                             ((alpha_b - porosity) * (1.0 - alpha_b) / K_S +
                              porosity * storage_compressibility + beta_FS);

        auto const b = (_process_data.initialisationLoadFactor(t) *
                    _process_data.specific_body_force).eval();

        // bodyforce-driven Darcy flow
        local_rhs.noalias() +=
            dNdx_p.transpose() * rho_fr * rho_fr * K_over_mu * b * w;

        // density dependence on pressure evaluated for Darcy-term,
        // for laplace and storage terms this dependence is neglected (as is
        // done for monolithic too)
        add_p_derivative.noalias() += rho_fr * beta_p * dNdx_p.transpose() *
                                      K_over_mu *
                                      (dNdx_p * p - 2.0 * rho_fr * b) * N_p * w;

        if (!fixed_stress_over_time_step)
        {
            auto const& eps = _ip_data[ip].eps;
            auto const& eps_prev = _ip_data[ip].eps_prev;
            const double eps_v_dot =
                (Invariants::trace(eps) - Invariants::trace(eps_prev)) / dt;

            // Constant portion of strain rate term:
            double const strain_rate_b =
                alpha_b * eps_v_dot -
                beta_FS * _ip_data[ip].strain_rate_variable;

            local_rhs.noalias() -= strain_rate_b * rho_fr * N_p * w;
        }
        else
        {
            // Constant portion of strain rate term:
            local_rhs.noalias() -=
                alpha_b * _ip_data[ip].strain_rate_variable * rho_fr * N_p * w;
        }
    }
    local_Jac.noalias() = laplace + storage / dt + add_p_derivative;

    local_rhs.noalias() -= laplace * p + storage * (p - p_prev) / dt;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    auto const p = local_x.template segment<pressure_size>(pressure_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, displacement_size>>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<displacement_size>>(
            local_b_data, displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& medium = _process_data.media_map.getMedium(_element.getID());
    auto const& solid = medium->phase(MaterialPropertyLib::PhaseName::Solid);
    auto const& fluid = fluidPhase(*medium);
    MPL::VariableArray vars;
    auto& phase_pressure = medium->hasPhase(MaterialPropertyLib::PhaseName::Gas)
                               ? vars.gas_phase_pressure
                               : vars.liquid_phase_pressure;

    auto const T_ref =
        medium->property(MPL::PropertyType::reference_temperature)
            .template value<double>(vars, x_position, t, dt);
    vars.temperature = T_ref;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N_p = _ip_data[ip].N_p;

        x_position = {
            std::nullopt, _element.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, N_u))};

        auto const x_coord = x_position.getCoordinates().value()[0];
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;

        phase_pressure = N_p.dot(p);

        auto const alpha = usesTotalStress() ? 0. : medium->property(MPL::PropertyType::biot_coefficient)
                               .template value<double>(vars, x_position, t, dt);
        auto const rho_sr =
            solid.property(MPL::PropertyType::density)
                .template value<double>(vars, x_position, t, dt);
        auto const porosity =
            medium->property(MPL::PropertyType::porosity)
                .template value<double>(vars, x_position, t, dt);

        auto const rho_fr = MaterialPropertyLib::getFluidDensity(
            t, dt, x_position, fluid, vars);

        auto const b = (_process_data.initialisationLoadFactor(t) *
                    _process_data.specific_body_force).eval();
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::kelvin_vector_dimensions(
                DisplacementDim)>::identity2;

        eps.noalias() = B * u;
        eps -= _ip_data[ip].birth_strain;
        vars.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        auto C = _ip_data[ip].updateConstitutiveRelation(vars, t, x_position,
                                                         dt, u, T_ref);

        local_Jac.noalias() += B.transpose() * C * B * w;

        double p_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(p, N_p, p_at_xi);

        double const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
        local_rhs.noalias() -=
            (B.transpose() * (sigma_eff - alpha * identity2 * p_at_xi) -
             N_u_op(N_u).transpose() * rho * b) *
            w;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(const double t, double const dt,
                                           Eigen::VectorXd const& local_x,
                                           Eigen::VectorXd const& local_x_prev,
                                           int const process_id,
                                           std::vector<double>& local_b_data,
                                           std::vector<double>& local_Jac_data)
{
    // For the equations with pressure
    if (process_id == _process_data.hydraulic_process_id)
    {
        assembleWithJacobianForPressureEquations(t, dt, local_x, local_x_prev,
                                                 local_b_data, local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(t, dt, local_x, local_b_data,
                                                local_Jac_data);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, DisplacementDim>::
    setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                 double const t,
                                 int const /*process_id*/)
{
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    auto const& medium = _process_data.media_map.getMedium(_element.getID());

    auto const p = local_x.template segment<pressure_size>(pressure_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto const& identity2 = Invariants::identity2;
    constexpr double dt = std::numeric_limits<double>::quiet_NaN();

    MPL::VariableArray vars;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        x_position = {
            std::nullopt, _element.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, N_u))};

        auto const x_coord = x_position.getCoordinates().value()[0];
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;
        eps -= _ip_data[ip].birth_strain;
        vars.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        if (_process_data.initial_stress.isTotalStress())
        {
            auto const& N_p = _ip_data[ip].N_p;
            auto const alpha_b =
                medium->property(MPL::PropertyType::biot_coefficient)
                    .template value<double>(vars, x_position, t, dt);

            auto& sigma_eff = _ip_data[ip].sigma_eff;
            sigma_eff.noalias() += alpha_b * N_p.dot(p) * identity2;
            _ip_data[ip].sigma_eff_prev.noalias() = sigma_eff;
        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, DisplacementDim>::
    postNonLinearSolverConcrete(Eigen::VectorXd const& local_x,
                                Eigen::VectorXd const& local_x_prev,
                                double const t, double const dt,
                                int const process_id)
{
    // Note: local_x and local_x_prev only contain the solutions of current
    // process in the staggered scheme. This has to be changed according to the
    // same two arguments in postTimestepConcrete.

    int const n_integration_points = _integration_method.getNumberOfPoints();

    auto const staggered_scheme_ptr =
        std::get_if<Staggered>(&_process_data.coupling_scheme);

    if (staggered_scheme_ptr &&
        process_id == _process_data.hydraulic_process_id)
    {
        if (!staggered_scheme_ptr->fixed_stress_over_time_step)
        {
            auto const p =
                local_x.template segment<pressure_size>(pressure_index);
            auto const p_prev =
                local_x_prev.template segment<pressure_size>(pressure_index);

            for (int ip = 0; ip < n_integration_points; ip++)
            {
                auto& ip_data = _ip_data[ip];

                auto const& N_p = ip_data.N_p;

                ip_data.strain_rate_variable = N_p.dot(p - p_prev) / dt;
            }
        }
    }

    if (!staggered_scheme_ptr ||
        process_id == _process_data.mechanics_related_process_id)
    {
        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        auto const& medium =
            _process_data.media_map.getMedium(_element.getID());

        auto const T_ref =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(MPL::EmptyVariableArray, x_position, t,
                                        dt);

        auto const u =
            local_x.template segment<displacement_size>(displacement_index);

        MPL::VariableArray vars;
        vars.temperature = T_ref;

        for (int ip = 0; ip < n_integration_points; ip++)
        {
            auto const& N_u = _ip_data[ip].N_u;
            auto const& dNdx_u = _ip_data[ip].dNdx_u;

            x_position = {std::nullopt, _element.getID(),
                          MathLib::Point3d(NumLib::interpolateCoordinates<
                                           ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                              _element, N_u))};

            auto const x_coord = x_position.getCoordinates().value()[0];
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunctionDisplacement::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx_u, N_u, x_coord,
                                                     _is_axially_symmetric);

            auto& eps = _ip_data[ip].eps;
            eps.noalias() = B * u;
        eps -= _ip_data[ip].birth_strain;
            vars.mechanical_strain.emplace<
                MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(eps);

            _ip_data[ip].updateConstitutiveRelation(vars, t, x_position, dt, u,
                                                    T_ref);
            if (!usesTotalStress() && _process_data.mechanical_drainage && _process_data.mechanical_drainage->active(t))
            {
                auto const p = local_x.template segment<pressure_size>(pressure_index);
                double const alpha = medium->property(MPL::PropertyType::biot_coefficient)
                    .template value<double>(vars, x_position, t, dt);
                double const porosity = medium->property(MPL::PropertyType::porosity)
                    .template value<double>(vars, x_position, t, dt);
                auto const C_el = _ip_data[ip].computeElasticTangentStiffness(t, x_position, dt, T_ref);
                double const K_d = _ip_data[ip].solid_material.getBulkModulus(t, x_position, &C_el);
                auto const material = _process_data.mechanical_drainage->materialAt(getMaterialID(),t)
                    .resolved(alpha, porosity, K_d);
                auto const response = _ip_data[ip].mechanical_drainage.evaluate(
                    _ip_data[ip].N_p.dot(p), Invariants::trace(eps), alpha, material);
                _ip_data[ip].mechanical_pressure = response.pressure;
                _ip_data[ip].mechanical_excess = response.excess;
                _ip_data[ip].mechanical_fluid_content = response.fluid_content;
            }

        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::postTimestepConcrete(Eigen::VectorXd const& local_x,
                                           Eigen::VectorXd const& local_x_prev,
                                           double const t, double const dt,
                                           int const process_id)
{
    auto const staggered_scheme_ptr =
        std::get_if<Staggered>(&_process_data.coupling_scheme);

    if (staggered_scheme_ptr &&
        process_id == _process_data.hydraulic_process_id)
    {
        if (staggered_scheme_ptr->fixed_stress_over_time_step)
        {
            auto const fixed_stress_stabilization_parameter =
                staggered_scheme_ptr->fixed_stress_stabilization_parameter;

            auto const p =
                local_x.template segment<pressure_size>(pressure_index);
            auto const p_prev =
                local_x_prev.template segment<pressure_size>(pressure_index);

            ParameterLib::SpatialPosition x_position;
            x_position.setElementID(_element.getID());

            auto const& solid_material =
                MaterialLib::Solids::selectSolidConstitutiveRelation(
                    _process_data.solid_materials, _process_data.material_ids,
                    _element.getID());

            auto const& medium =
                _process_data.media_map.getMedium(_element.getID());
            MPL::VariableArray vars;

            auto const T_ref =
                medium->property(MPL::PropertyType::reference_temperature)
                    .template value<double>(vars, x_position, t, dt);
            vars.temperature = T_ref;

            int const n_integration_points =
                _integration_method.getNumberOfPoints();
            for (int ip = 0; ip < n_integration_points; ip++)
            {
                auto& ip_data = _ip_data[ip];

                auto const& N_p = ip_data.N_p;

                x_position = {
                    std::nullopt, _element.getID(),
                    MathLib::Point3d(
                        NumLib::interpolateCoordinates<
                            ShapeFunctionPressure, ShapeMatricesTypePressure>(
                            _element, N_p))};

                auto const& eps = ip_data.eps;
                auto const& eps_prev = ip_data.eps_prev;
                const double eps_v_dot =
                    (Invariants::trace(eps) - Invariants::trace(eps_prev)) / dt;

                auto const C_el = ip_data.computeElasticTangentStiffness(
                    t, x_position, dt, T_ref);
                auto const K_S =
                    solid_material.getBulkModulus(t, x_position, &C_el);

                auto const alpha_b =
                    medium->property(MPL::PropertyType::biot_coefficient)
                        .template value<double>(vars, x_position, t, dt);

                ip_data.strain_rate_variable =
                    eps_v_dot - fixed_stress_stabilization_parameter * alpha_b *
                                    N_p.dot(p - p_prev) / dt / K_S;
            }
        }
    }

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data[ip].pushBackState();
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::size_t HydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::setIPDataInitialConditions(std::string_view const name,
                                                 double const* values,
                                                 int const integration_order)
{
    if (integration_order !=
        static_cast<int>(_integration_method.getIntegrationOrder()))
    {
        OGS_FATAL(
            "Setting integration point initial conditions; The integration "
            "order of the local assembler for element {:d} is different from "
            "the integration order in the initial condition.",
            _element.getID());
    }

    if (name == "sigma")
    {
        if (_process_data.initial_stress.value != nullptr)
        {
            OGS_FATAL(
                "Setting initial conditions for stress from integration "
                "point data and from a parameter '{:s}' is not possible "
                "simultaneously.",
                _process_data.initial_stress.value->name);
        }

        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::sigma_eff);
    }

    if (name == "epsilon")
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::eps);
    }

    if (name == "strain_rate_variable")
    {
        return ProcessLib::setIntegrationPointScalarData(
            values, _ip_data, &IpData::strain_rate_variable);
    }

    return 0;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double>
HydroMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                             DisplacementDim>::getSigma() const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IpData::sigma_eff);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double>
HydroMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                             DisplacementDim>::getEpsilon() const
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<double> ip_epsilon_values;
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
        ip_epsilon_values, n_integration_points, kelvin_vector_size);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& eps = _ip_data[ip].eps;
        cache_mat.row(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps);
    }

    return ip_epsilon_values;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
std::vector<double>
HydroMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                             DisplacementDim>::getStrainRateVariable() const
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<double> ip_strain_rate_variables(n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        ip_strain_rate_variables[ip] = _ip_data[ip].strain_rate_variable;
    }

    return ip_strain_rate_variables;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                  ShapeFunctionPressure, DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& /*local_x_prev*/)
{
    auto const p = local_x.template segment<pressure_size>(pressure_index);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, p,
                         *_process_data.pressure_interpolated);

    int const elem_id = _element.getID();

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto const& medium = _process_data.media_map.getMedium(elem_id);
    MPL::VariableArray vars;
    auto& phase_pressure = medium->hasPhase(MaterialPropertyLib::PhaseName::Gas)
                               ? vars.gas_phase_pressure
                               : vars.liquid_phase_pressure;

    SymmetricTensor k_sum = SymmetricTensor::Zero(KelvinVectorSize);
    auto sigma_eff_sum = MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
        Eigen::Matrix<double, 3, 3>::Zero());

    auto const& identity2 = Invariants::identity2;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& eps = _ip_data[ip].eps;
        sigma_eff_sum += _ip_data[ip].sigma_eff;

        ParameterLib::SpatialPosition x_position = {
            std::nullopt, _element.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, _ip_data[ip].N_u))};

        auto const alpha = usesTotalStress() ? 0. : medium->property(MPL::PropertyType::biot_coefficient)
                               .template value<double>(vars, x_position, t, dt);
        double const p_int_pt = _ip_data[ip].N_p.dot(p);

        phase_pressure = p_int_pt;

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            auto const sigma_total =
                (_ip_data[ip].sigma_eff - alpha *
                 (_process_data.use_positive_pore_pressure_in_mechanics
                      ? std::max(p_int_pt, 0.) : p_int_pt) * identity2).eval();
            vars.total_stress.emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    sigma_total));
        }
        // For strain dependent permeability
        vars.volumetric_strain = Invariants::trace(eps);
        vars.equivalent_plastic_strain =
            _ip_data[ip].material_state_variables->getEquivalentPlasticStrain();
        vars.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        k_sum += MPL::getSymmetricTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(vars, x_position, t, dt));
    }

    Eigen::Map<Eigen::VectorXd>(
        &(*_process_data.permeability)[elem_id * KelvinVectorSize],
        KelvinVectorSize) = k_sum / n_integration_points;

    Eigen::Matrix<double, 3, 3, 0, 3, 3> const sigma_avg =
        MathLib::KelvinVector::kelvinVectorToTensor(sigma_eff_sum) /
        n_integration_points;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(sigma_avg);

    Eigen::Map<Eigen::Vector3d>(
        &(*_process_data.principal_stress_values)[elem_id * 3], 3) =
        e_s.eigenvalues();

    auto eigen_vectors = e_s.eigenvectors();

    for (auto i = 0; i < 3; i++)
    {
        Eigen::Map<Eigen::Vector3d>(
            &(*_process_data.principal_stress_vector[i])[elem_id * 3], 3) =
            eigen_vectors.col(i);
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
unsigned HydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::getNumberOfIntegrationPoints() const
{
    return _integration_method.getNumberOfPoints();
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
int HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                 ShapeFunctionPressure,
                                 DisplacementDim>::getMaterialID() const
{
    return _process_data.material_ids == nullptr
               ? 0
               : (*_process_data.material_ids)[_element.getID()];
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
typename MaterialLib::Solids::MechanicsBase<
    DisplacementDim>::MaterialStateVariables const&
HydroMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
                             DisplacementDim>::
    getMaterialStateVariablesAt(unsigned integration_point) const
{
    return *_ip_data[integration_point].material_state_variables;
}
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::acceptMechanicalDrainage(std::size_t const /*mesh_item_id*/, double const t, GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table)
{
    auto const& settings = _process_data.mechanical_drainage;
    if (!settings || (settings->active(t) && !settings->isSafetyReference(t))) { return; }
    if (settings->isSafetyReference(t) && settings->materials.at(_process_data.material_ids ? (*_process_data.material_ids)[_element.getID()] : 0).usesTotalStress()) { return; }
    auto const local = x.get(NumLib::getIndices(_element.getID(), dof_table));
    auto const p = Eigen::Map<Eigen::VectorXd const>(local.data() + pressure_index, pressure_size);
    auto const u = Eigen::Map<Eigen::VectorXd const>(local.data() + displacement_index, displacement_size);
    for (auto& ip : _ip_data)
    {
        ParameterLib::SpatialPosition const pos{std::nullopt, _element.getID(),
            MathLib::Point3d(NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                ShapeMatricesTypeDisplacement>(_element, ip.N_u))};
        auto const B = LinearBMatrix::computeBMatrix<DisplacementDim,
            ShapeFunctionDisplacement::NPOINTS, typename BMatricesType::BMatrixType>(
                ip.dNdx_u, ip.N_u, pos.getCoordinates().value()[0], _is_axially_symmetric);
        double const primary=ip.N_p.dot(p);
        if (settings->isSafetyReference(t))
        {
            // Before the first strict stage, ordinary positive-Biot output may
            // be clipped. Capture the signed primary fluid state, never that
            // clipped diagnostic. A prior strict stage already has actual p.
            if (settings->reference_time == settings->safety_reference_time)
            { ip.mechanical_pressure = primary; }
            ip.mechanical_drainage.rebaseAccepted(Invariants::trace((B*u).eval()),primary,ip.mechanical_pressure);
            ip.mechanical_fluid_content=0.;
        }
        else
        {
            ip.mechanical_pressure = primary;
            ip.mechanical_drainage.accept(t, settings->reference_time,
                Invariants::trace((B*u).eval()), ip.mechanical_pressure);
        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure, int DisplacementDim>
void HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,ShapeFunctionPressure,DisplacementDim>::
preTimestepConcrete(std::vector<double> const& local_x,double const t,double const dt)
{
    for (auto const& birth:_process_data.construction_births)
    {
        if (birth.material_id!=getMaterialID() || !(NumLib::Time{t}>NumLib::Time{birth.time}) ||
            NumLib::Time{t-dt}>NumLib::Time{birth.time} || applied_construction_births_.contains(birth.time)) { continue; }
        if (birth.reset)
        {
            auto const u=Eigen::Map<Eigen::VectorXd const>(local_x.data()+displacement_index,displacement_size);
            for (auto& ip:_ip_data)
            {
                ParameterLib::SpatialPosition const pos{std::nullopt,_element.getID(),MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunctionDisplacement,ShapeMatricesTypeDisplacement>(_element,ip.N_u))};
                auto const B=LinearBMatrix::computeBMatrix<DisplacementDim,ShapeFunctionDisplacement::NPOINTS,
                    typename BMatricesType::BMatrixType>(ip.dNdx_u,ip.N_u,pos.getCoordinates().value()[0],_is_axially_symmetric);
                ip.birth_strain=B*u;ip.eps.setZero();ip.eps_prev.setZero();ip.sigma_eff.setZero();ip.sigma_eff_prev.setZero();
                ip.material_state_variables=ip.solid_material.createMaterialStateVariables();
                ip.solid_material.initializeInternalStateVariables(birth.time,pos,*ip.material_state_variables);
                ip.pushBackState();
            }
        }
        applied_construction_births_.insert(birth.time);
    }
}

}  // namespace HydroMechanics
}  // namespace ProcessLib
