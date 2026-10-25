// SPDX-License-Identifier: BSD-3-Clause
#include "BacktrackingNewtonStrategy.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "DampingPolicy.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NonlinearSystem.h"
#include "NumLib/Exceptions.h"

namespace NumLib
{
BacktrackingNewtonStrategy::BacktrackingNewtonStrategy(
    double initial,double contraction,double minimum,double armijo,
    unsigned max_trials,double merit_tolerance, bool residual_screening)
    :initial_(initial),contraction_(contraction),minimum_(minimum),armijo_(armijo),
     merit_tolerance_(merit_tolerance),max_trials_(max_trials),
     residual_screening_(residual_screening)
{
    if (!std::isfinite(initial) || !std::isfinite(contraction) ||
        !std::isfinite(minimum) || !std::isfinite(armijo) ||
        !std::isfinite(merit_tolerance) || !(0.<minimum && minimum<=initial && initial<=1.) ||
        !(0.<contraction && contraction<1.) || !(0.<armijo && armijo<1.) ||
        max_trials==0 || !(0.<merit_tolerance))
    { OGS_FATAL("Invalid Newton Backtracking settings."); }
}
StepResult BacktrackingNewtonStrategy::applyStep(
    GlobalVector const& x,GlobalVector const& direction,
    GlobalVector const&,GlobalMatrix const& jacobian,GlobalVector& next,
    NewtonStepContext& context,int)
{
#if defined(USE_PETSC) || defined(USE_LIS)
    OGS_FATAL("Backtracking currently requires serial Eigen monolithic solves.");
#else
    namespace LA=MathLib::LinAlg;
    if (!context.current_x || context.current_x->size()!=1 || context.process_id!=0)
    { OGS_FATAL("Backtracking currently requires one monolithic equation system."); }
    if (residual_screening_ && !context.sys.supportsResidualScreening())
    { OGS_FATAL("ResidualThenJacobian is unsupported for this process/configuration."); }
    // Recover the physical Jacobian/residual from the equation system rather
    // than the arguments: Eigen equilibration modifies the solve buffers.
    GlobalMatrix trial_jacobian{jacobian};
    GlobalVector trial_residual{x},scratch{x},base{x};
    context.sys.computeKnownSolutions(base,0);
    context.sys.applyKnownSolutions(base);
    Eigen::VectorXd row_scales;
    unsigned assembly_calls = 0;
    double assembly_seconds = 0.;
    unsigned residual_calls = 0;
    unsigned validated_candidates = 0;
    auto evaluate=[&](GlobalVector& trial,bool const establish_scales,
                          bool const reuse_assembly=false,
                          bool const residual_only=false) -> std::optional<double>
    {
        auto current=*context.current_x;current[0]=&trial;
        context.sys.computeKnownSolutions(trial,0);
        context.sys.applyKnownSolutions(trial);
        if (!reuse_assembly)
        {
            BaseLib::RunTime assembly_timer;
            assembly_timer.start();
            ++assembly_calls;
            try
            {
                if (residual_only)
                {
                    ++residual_calls;
                    context.sys.assembleResidual(current,context.x_prev,0);
                }
                else { context.sys.assemble(current,context.x_prev,0); }
            }
            catch (AssemblyException const& error)
            {
                assembly_seconds += assembly_timer.elapsed();
                WARN("Backtracking rejected an inadmissible assembly: {}",error.what());
                return std::nullopt;
            }
            assembly_seconds += assembly_timer.elapsed();
        }
        context.sys.getResidual(trial,*context.x_prev[0],trial_residual);
        context.sys.getJacobian(trial_jacobian);
        if (context.initial_residual) { LA::axpy(trial_residual,-1.,*context.initial_residual); }
        scratch.setZero();
        context.sys.applyKnownSolutionsNewton(trial_jacobian,trial_residual,trial,scratch);
        if (establish_scales)
        {
            auto const& matrix=trial_jacobian.getRawMatrix();
            row_scales=Eigen::VectorXd::Zero(matrix.rows());
            for (Eigen::Index k=0;k<matrix.outerSize();++k)
            {
                for (GlobalMatrix::RawMatrixType::InnerIterator it(matrix,k);it;++it)
                {
                    double const characteristic=std::max(1.,std::abs(base[it.col()]));
                    row_scales[it.row()]=std::max(row_scales[it.row()],std::abs(it.value())*characteristic);
                }
            }
            row_scales=row_scales.cwiseMax(1e-30);
        }
        auto const scaled=(trial_residual.getRawVector().array()/row_scales.array()).matrix();
        double const merit=.5*scaled.squaredNorm();
        return std::isfinite(merit)?std::optional<double>{merit}:std::nullopt;
    };
    // Newton has already assembled at x. Its physical buffers remain in
    // the equation system even when the linear solver equilibrates its copies.
    // Reuse only if applying the current essential conditions changed no DOF;
    // otherwise the constrained base is a different state and needs assembly.
    bool const reuse_base =
        (base.getRawVector().array() == x.getRawVector().array()).all();
    auto const initial_merit=evaluate(base,true,reuse_base);
    DBUG("Backtracking base assembly reused: {}", reuse_base);
    if (!initial_merit) { return {false,0.,true}; }
    double step=initial_;
    if (_damping_policy) { step=_damping_policy->apply(direction,x,step); }
    if (!std::isfinite(step) || step<minimum_) { return {false,step,true}; }
    unsigned trials=0;
    auto const accepted=backtrackArmijo(step,contraction_,minimum_,armijo_,max_trials_,
        *initial_merit,.5*merit_tolerance_*merit_tolerance_,[&](double const alpha)
        {
            ++trials;LA::copy(base,next);LA::axpy(next,-alpha,direction);
            auto merit = evaluate(next,false,false,residual_screening_);
            if (residual_screening_ && merit &&
                (*merit <= (1.-armijo_*alpha)*(*initial_merit) ||
                 *merit <= .5*merit_tolerance_*merit_tolerance_))
            {
                // Screening cannot certify a tangent or the complete fixed-N
                // substepped update. Fully re-evaluate before any acceptance.
                ++validated_candidates;
                merit = evaluate(next,false);
            }
            return merit;
        });
    INFO("[time] Backtracking assemblies: {}, total {:g} s; base reused: {}.",
         assembly_calls, assembly_seconds, reuse_base);
    if (residual_screening_)
    {
        INFO("Residual screening: {} residual trials, {} full candidate validations.",
             residual_calls, validated_candidates);
    }
    if (!accepted)
    {
        LA::copy(x,next);
        WARN("Backtracking exhausted {} trials; reject the global increment.",trials);
        return {false,0.,true};
    }
    INFO("Backtracking accepted alpha={} after {} trials; scaled merit {} -> {}.",
        accepted->first,trials,*initial_merit,accepted->second);
    return {true,accepted->first,true};
#endif
}
} // namespace NumLib
