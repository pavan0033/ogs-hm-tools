// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <set>

#include <memory>
#include <vector>

#include "HydroMechanicsProcessData.h"
#ifdef OGS_USE_MFRONT
#include "IntegrationPointCheckpoint.h"
#include "MaterialLib/SolidModels/MFront/MFront.h"
#endif
#include "LocalAssemblerInterface.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"

namespace ProcessLib
{
namespace HydroMechanics
{
namespace MPL = MaterialPropertyLib;

template <typename BMatricesType, typename ShapeMatricesTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;
    typename BMatricesType::KelvinVectorType birth_strain = BMatricesType::KelvinVectorType::Zero();

    typename ShapeMatricesTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatricesTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;
    double integration_weight = 0;

    // previous pressure rate for the fixed stress splitting
    // approach in the staggered scheme.
    // TODO: disable in monolithic scheme to save memory.
    double strain_rate_variable = 0.0;
    MechanicalDrainageState mechanical_drainage;
    double mechanical_pressure = 0.;
    double mechanical_excess = 0.;
    double mechanical_fluid_content = 0.;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }

    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>
    computeElasticTangentStiffness(
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        double const temperature)
    {
#ifdef OGS_USE_MFRONT
        if (auto const* mfront = dynamic_cast<MaterialLib::Solids::MFront::MFront<DisplacementDim> const*>(&solid_material))
        {
            if (auto const elastic = mfront->getMCElasticStiffness(t, x_position))
            {
                return *elastic;
            }
            if (auto const elastic = mfront->getHSElasticStiffness(t, x_position, sigma_eff))
            {
                return *elastic;
            }
        }
#endif
        namespace MPL = MaterialPropertyLib;

        MPL::VariableArray variable_array;
        MPL::VariableArray variable_array_prev;

        auto const null_state = solid_material.createMaterialStateVariables();
        solid_material.initializeInternalStateVariables(t, x_position,
                                                        *null_state);

        using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

        variable_array.stress.emplace<KV>(KV::Zero());
        variable_array.mechanical_strain.emplace<KV>(KV::Zero());
        variable_array.temperature = temperature;

        variable_array_prev.stress.emplace<KV>(KV::Zero());
        variable_array_prev.mechanical_strain.emplace<KV>(KV::Zero());
        variable_array_prev.temperature = temperature;

        auto&& solution =
            solid_material.integrateStress(variable_array_prev, variable_array,
                                           t, x_position, dt, *null_state);

        if (!solution)
        {
            OGS_FATAL("Computation of elastic tangent stiffness failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C =
            std::move(std::get<2>(*solution));

        return C;
    }

    template <typename DisplacementVectorType>
    typename BMatricesType::KelvinMatrixType updateConstitutiveRelation(
        MPL::VariableArray const& variable_array,
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        DisplacementVectorType const& /*u*/,
        double const T, bool const compute_tangent = true)
    {
        MaterialPropertyLib::VariableArray variable_array_prev;
        variable_array_prev.stress
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_eff_prev);
        variable_array_prev.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_prev);
        variable_array_prev.temperature = T;

        auto solution = [&]()
        {
#ifdef OGS_USE_MFRONT
            if (!compute_tangent)
            {
                if (auto const* mfront = dynamic_cast<MaterialLib::Solids::MFront::MFront<DisplacementDim> const*>(&solid_material))
                {
                    return mfront->integrateStressForResidual(variable_array_prev,
                        variable_array, t, x_position, dt, *material_state_variables);
                }
            }
#endif
            return solid_material.integrateStress(variable_array_prev,
                variable_array, t, x_position, dt, *material_state_variables);
        }();

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables, C) = std::move(*solution);

        return C;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used for the extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N_u;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
class HydroMechanicsLocalAssembler
    : public LocalAssemblerInterface<DisplacementDim>
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    using SymmetricTensor = Eigen::Matrix<double, KelvinVectorSize, 1>;

    static constexpr auto& N_u_op = MathLib::eigenBlockMatrixView<
        DisplacementDim,
        typename ShapeMatricesTypeDisplacement::NodalRowVectorType>;

    HydroMechanicsLocalAssembler(HydroMechanicsLocalAssembler const&) = delete;
    HydroMechanicsLocalAssembler(HydroMechanicsLocalAssembler&&) = delete;

    HydroMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        HydroMechanicsProcessData<DisplacementDim>& process_data);

    std::vector<std::string> captureIntegrationPointCheckpoint() const override
    {
#ifdef OGS_USE_MFRONT
        std::vector<std::string> result;
        result.reserve(_ip_data.size());
        for (auto const& ip : _ip_data)
        { result.push_back(captureHMIntegrationPoint<DisplacementDim>(ip, usesTotalStress())); }
        return result;
#else
        OGS_FATAL("HM integration-point checkpoint requires MFront.");
#endif
    }
    void validateIntegrationPointCheckpoint(std::vector<std::string> const& capsules) const override
    {
#ifdef OGS_USE_MFRONT
        if (capsules.size()!=_ip_data.size()) { OGS_FATAL("HM checkpoint IP count mismatch."); }
        for (std::size_t i=0;i<capsules.size();++i)
        { (void)decodeHMIntegrationPoint<DisplacementDim>(_ip_data[i],capsules[i],usesTotalStress()); }
#else
        OGS_FATAL("HM checkpoint requires MFront.");
#endif
    }
    void restoreIntegrationPointCheckpoint(std::vector<std::string> const& capsules) override
    {
#ifdef OGS_USE_MFRONT
        if (capsules.size() != _ip_data.size())
        { OGS_FATAL("HM checkpoint integration-point count mismatch."); }
        // Validate the whole element before mutating any receiving IP.
        std::vector<DecodedHMIntegrationPoint<DisplacementDim>> decoded;
        decoded.reserve(capsules.size());
        for (std::size_t i = 0; i < capsules.size(); ++i)
        { decoded.push_back(decodeHMIntegrationPoint<DisplacementDim>(_ip_data[i], capsules[i], usesTotalStress())); }
        for (std::size_t i = 0; i < capsules.size(); ++i)
        { commitHMIntegrationPoint<DisplacementDim>(_ip_data[i], std::move(decoded[i])); }
#else
        OGS_FATAL("HM integration-point checkpoint requires MFront.");
#endif
    }

    /// Returns number of read integration points.
    double getInitialisationYieldResidual(double const factor=1.) const override
    {
        if constexpr (DisplacementDim != 2)
        { OGS_FATAL("MC initialisation gate requires plane strain."); }
        else
        {
            auto const& settings = *_process_data.initialisation_gate;
            int const id = _process_data.material_ids ? (*_process_data.material_ids)[_element.getID()] : 0;
            auto const base_material=settings.materials.at(id);
            double worst = -std::numeric_limits<double>::infinity();
            for (auto const& ip : _ip_data)
            {
                auto material = base_material;
                if (auto const parameter=settings.cohesion_parameters.find(id); parameter!=settings.cohesion_parameters.end())
                {
                    ParameterLib::SpatialPosition const pos{std::nullopt, _element.getID(),
                        MathLib::Point3d(NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                            ShapeMatricesTypeDisplacement>(_element, ip.N_u))};
                    material.cohesion=(*parameter->second)(0.,pos)[0];
                }
                double const cap=material.friction_angle>0. ? std::min(100.,.1*material.cohesion/std::tan(material.friction_angle*std::numbers::pi/180.)) : 100.;
                if (!std::isfinite(material.cohesion) || material.cohesion<=0. || material.cutoff>cap)
                { OGS_FATAL("Spatial MC cohesion/cutoff is inadmissible at an integration point."); }
                if (factor!=1.) { material.cohesion/=factor;material.friction_angle=std::atan(std::tan(material.friction_angle*std::numbers::pi/180.)/factor)*180./std::numbers::pi; }
                if (!ip.eps.allFinite() || !ip.sigma_eff.allFinite() ||
                    !std::isfinite(ip.material_state_variables->getEquivalentPlasticStrain()))
                { OGS_FATAL("Initialisation gate: nonfinite constitutive state."); }
                worst = std::max(worst, mcAdmissibilityResidual(ip.sigma_eff, material));
            }
            return worst;
        }
    }

    std::size_t setIPDataInitialConditions(
        std::string_view const name,
        double const* values,
        int const integration_order) override;

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_x_prev*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "HydroMechanicsLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override;

    void assembleWithJacobianForStaggeredScheme(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, int const process_id,
        std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override;

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = _ip_data[ip];

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, _element.getID(),
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunctionPressure,
                                                   ShapeMatricesTypePressure>(
                        _element, ip_data.N_p))};

            /// Set initial stress from parameter.
            if (_process_data.initial_stress.value)
            {
                ip_data.sigma_eff =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector<
                        DisplacementDim>((*_process_data.initial_stress.value)(
                        std::numeric_limits<
                            double>::quiet_NaN() /* time independent */,
                        x_position));
            }

            double const t = 0;  // TODO (naumov) pass t from top
            ip_data.solid_material.initializeInternalStateVariables(
                t, x_position, *ip_data.material_state_variables);

            ip_data.pushBackState();
        }
    }

    void preTimestepConcrete(std::vector<double> const& local_x,
                             double const t,double const dt) override;

    void postTimestepConcrete(Eigen::VectorXd const& local_x,
                              Eigen::VectorXd const& local_x_prev,
                              double const t, double const dt,
                              int const process_id) override;

    void computeSecondaryVariableConcrete(
        double const t, double const dt, Eigen::VectorXd const& local_xs,
        Eigen::VectorXd const& local_x_prev) override;

    void postNonLinearSolverConcrete(Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_prev,
                                     double const t, double const dt,
                                     int const process_id) override;

    void setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                      double const t,
                                      int const process_id) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_u = _secondary_data.N_u[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_u.data(), N_u.size());
    }

    // TODO (naumov) This method is same as getIntPtSigma but for arguments and
    // the ordering of the cache_mat.
    // There should be only one.
    void acceptMechanicalDrainage(std::size_t mesh_item_id, double t, GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table) override;
    bool usesTotalStress() const
    {
        return _process_data.mechanical_drainage &&
            _process_data.mechanical_drainage->materials.at(getMaterialID()).usesTotalStress();
    }
    std::vector<double> getMechanicalStressBasis() const override
    {
        return std::vector<double>(_ip_data.size(), usesTotalStress() ? 1. : 0.);
    }
    std::vector<double> getMechanicalPressureValid() const override
    {
        return std::vector<double>(_ip_data.size(), usesTotalStress() ? 0. : 1.);
    }
    std::vector<double> getMechanicalPressure() const override
    {
        if (usesTotalStress()) { return std::vector<double>(_ip_data.size(), std::numeric_limits<double>::quiet_NaN()); }
        std::vector<double> result;
        for (auto const& ip : _ip_data) { result.push_back(ip.mechanical_pressure); }
        return result;
    }
    std::vector<double> getAppliedBiotPressure() const override
    {
        auto result = getMechanicalPressure();
        if (_process_data.use_positive_pore_pressure_in_mechanics)
        {
            auto const& settings = _process_data.mechanical_drainage;
            for (std::size_t i = 0; i < result.size(); ++i)
            {
                auto const& ip = _ip_data[i];
                if (settings && ip.mechanical_drainage.hasReference())
                {
                    auto const& material = settings->materials.at(getMaterialID());
                    result[i] = material.appliedPressure(result[i], ip.mechanical_excess,
                        ip.mechanical_drainage.referencePressure(), true);
                }
                else { result[i] = std::max(result[i], 0.); }
            }
        }
        return result;
    }
    std::vector<double> getMechanicalExcess() const override
    {
        if (usesTotalStress()) { return std::vector<double>(_ip_data.size(), std::numeric_limits<double>::quiet_NaN()); }
        std::vector<double> result;
        for (auto const& ip : _ip_data) { result.push_back(ip.mechanical_excess); }
        return result;
    }
    std::vector<double> getMechanicalFluidContent() const override
    {
        if (usesTotalStress()) { return std::vector<double>(_ip_data.size(), std::numeric_limits<double>::quiet_NaN()); }
        std::vector<double> result;
        for (auto const& ip : _ip_data) { result.push_back(ip.mechanical_fluid_content); }
        return result;
    }
    std::vector<double> getMechanicalReferenceStrain() const override
    {
        std::vector<double> result;
        for (auto const& ip : _ip_data) { result.push_back(ip.mechanical_drainage.referenceStrain()); }
        return result;
    }
    std::vector<double> getMechanicalReferencePressure() const override
    {
        if (usesTotalStress()) { return std::vector<double>(_ip_data.size(), std::numeric_limits<double>::quiet_NaN()); }
        std::vector<double> result;
        for (auto const& ip : _ip_data) { result.push_back(ip.mechanical_drainage.referencePressure()); }
        return result;
    }
    std::vector<double> getMechanicalActualReferencePressure() const override
    {
        if (usesTotalStress()) { return std::vector<double>(_ip_data.size(),std::numeric_limits<double>::quiet_NaN()); }
        std::vector<double> values;
        for (auto const& ip:_ip_data) { values.push_back(ip.mechanical_drainage.referenceActualPressure()); }
        return values;
    }
    std::vector<double> getSigma() const override;

    std::vector<double> getEpsilon() const override;

    std::vector<double> getStrainRateVariable() const override;

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma_eff, cache);
    }

    std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::eps, cache);
    }

    int getNumberOfVectorElementsForDeformation() const override
    {
        return displacement_size;
    }

private:
    /**
     * Assemble local matrices and vectors arise from the linearized discretized
     * weak form of the residual of the momentum balance equation,
     *      \f[
     *            \nabla (\sigma - \alpha_b p \mathrm{I}) = f
     *      \f]
     * where \f$ \sigma\f$ is the effective stress tensor, \f$p\f$ is the pore
     * pressure, \f$\alpha_b\f$ is the Biot constant, \f$\mathrm{I}\f$ is the
     * identity tensor, and \f$f\f$ is the body force.
     *
     * @param t               Time
     * @param dt              Time increment
     * @param local_x         Nodal values of \f$x\f$ of an element of all
     * coupled processes.
     * @param local_b_data    Right hand side vector of an element.
     * @param local_Jac_data  Element Jacobian matrix for the Newton-Raphson
     *                        method.
     */
    void assembleWithJacobianForDeformationEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data);

    /**
     * Assemble local matrices and vectors arise from the linearized discretized
     * weak form of the residual of the mass balance equation of single phase
     * flow,
     *      \f[
     *          \alpha \cdot{p} - \nabla (K (\nabla p + \rho g \nabla z) +
     *          \alpha_b \nabla \cdot \dot{u}  = Q
     *      \f]
     * where \f$ alpha\f$ is a coefficient may depend on storage or the fluid
     * density change, \f$ \rho\f$ is the fluid density, \f$g\f$ is the
     * gravitational acceleration, \f$z\f$ is the vertical coordinate, \f$u\f$
     * is the displacement, and \f$Q\f$ is the source/sink term.
     *
     * @param t               Time
     * @param dt              Time increment
     * @param local_x         Nodal values of \f$x\f$ of an element  of all
     * coupled processes.
     * @param local_x_prev    Nodal values of \f$x^{t-1}\f$ of an element  of
     * all coupled processes.
     * @param local_b_data    Right hand side vector of an element.
     * @param local_Jac_data  Element Jacobian matrix for the Newton-Raphson
     *                        method.
     */
    void assembleWithJacobianForPressureEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data);

    unsigned getNumberOfIntegrationPoints() const override;

    std::set<double> applied_construction_births_;
    int getMaterialID() const override;

    std::vector<double> getMaterialStateVariableInternalState(
        typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
            InternalVariable::WriteAccess const& values,
        int const& components) const override
    {
        return ProcessLib::getIntegrationPointDataMaterialStateVariables(
            _ip_data, &IpData::material_state_variables, values, components);
    }


    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const override;

private:
    HydroMechanicsProcessData<DisplacementDim>& _process_data;

    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeDisplacement,
                             ShapeMatricesTypePressure, DisplacementDim,
                             ShapeFunctionDisplacement::NPOINTS>;
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    NumLib::GenericIntegrationMethod const& _integration_method;
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;
    SecondaryData<
        typename ShapeMatricesTypeDisplacement::ShapeMatrices::ShapeType>
        _secondary_data;

    static const int pressure_index = 0;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib

#include "HydroMechanicsFEM-impl.h"
