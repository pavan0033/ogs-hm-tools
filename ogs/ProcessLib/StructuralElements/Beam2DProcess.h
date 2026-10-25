// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include "ProcessLib/Process.h"
#include "Beam2DLocalAssembler.h"
#include "PlasticFrameLocalAssembler.h"
namespace ProcessLib::StructuralElements
{
class PROCESSLIB_EXPORT Beam2DProcess final : public Process
{
public:
    Beam2DProcess(std::string name,MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
        unsigned integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&& variables,
        Beam2DProcessData data);
    bool isLinear() const override { return data_.axial_capacity==0.; }
private:
    void initializeConcreteProcess(NumLib::LocalToGlobalIndexMap const&,
        MeshLib::Mesh const&,unsigned) override;
    void assembleConcreteProcess(double,double,
        std::vector<GlobalVector*> const&,std::vector<GlobalVector*> const&,
        int,GlobalMatrix&,GlobalMatrix&,GlobalVector&) override;
    void assembleWithJacobianConcreteProcess(double,double,
        std::vector<GlobalVector*> const&,std::vector<GlobalVector*> const&,
        int,GlobalVector&,GlobalMatrix&) override;
    void computeSecondaryVariableConcrete(double,double,
        std::vector<GlobalVector*> const&,GlobalVector const&,int) override;
    void postAcceptedTimestepConcreteProcess(
        std::vector<GlobalVector*> const&,double,int) override;
    std::vector<GlobalIndexType> indices(MeshLib::Element const&) const;
    Beam2DProcessData const data_;
    std::vector<Beam2DLocalAssembler> assemblers_;
    std::vector<PlasticFrameLocalAssembler> plastic_assemblers_;
    std::vector<std::array<SectionState,2>> accepted_sections_;
    MeshLib::PropertyVector<double>* section_forces_=nullptr;
    MeshLib::PropertyVector<double>* plastic_strains_=nullptr;
    MeshLib::PropertyVector<double>* forces_=nullptr;
};
}
