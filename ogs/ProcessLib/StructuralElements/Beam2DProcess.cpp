// SPDX-License-Identifier: BSD-3-Clause
#include "Beam2DProcess.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
namespace ProcessLib::StructuralElements
{
Beam2DProcess::Beam2DProcess(std::string name,MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&& variables,
    Beam2DProcessData const data)
    :Process(std::move(name),mesh,std::move(jacobian),parameters,integration_order,
             std::move(variables),SecondaryVariableCollection{}),data_(data)
{
#ifdef USE_PETSC
    OGS_FATAL("BEAM_MECHANICS first scope requires serial Eigen.");
#endif
}
void Beam2DProcess::initializeConcreteProcess(NumLib::LocalToGlobalIndexMap const&,
    MeshLib::Mesh const& mesh,unsigned)
{
    if(mesh.isAxiallySymmetric()) { OGS_FATAL("BEAM_MECHANICS is planar only."); }
    for(auto const* element:mesh.getElements())
    {
        if(element->getDimension()!=1 || element->getNumberOfNodes()!=2)
        { OGS_FATAL("BEAM_MECHANICS requires two-node line elements."); }
        auto const delta=(element->getNode(1)->asEigenVector3d()-element->getNode(0)->asEigenVector3d()).eval();
        if(std::abs(delta[2])>1e-12 || delta.head<2>().norm()<=0.)
        { OGS_FATAL("Beam members must be nonzero length and lie in the xy plane."); }
        assemblers_.emplace_back(delta.head<2>(),data_);
        if(data_.axial_capacity>0.)
        {
            plastic_assemblers_.emplace_back(delta.head<2>(),data_.EA,data_.EI,data_.axial_capacity,data_.moment_capacity);
            accepted_sections_.emplace_back();
        }
    }
    if(data_.axial_capacity>0.)
    {
        section_forces_=MeshLib::getOrCreateMeshProperty<double>(_mesh,"beam_section_forces",MeshLib::MeshItemType::Cell,4);
        plastic_strains_=MeshLib::getOrCreateMeshProperty<double>(_mesh,"beam_plastic_deformation",MeshLib::MeshItemType::Cell,4);
    }
    forces_=MeshLib::getOrCreateMeshProperty<double>(_mesh,"beam_end_forces",MeshLib::MeshItemType::Cell,6);
}
std::vector<GlobalIndexType> Beam2DProcess::indices(MeshLib::Element const& element) const
{
    std::vector<GlobalIndexType> result;
    for(unsigned i=0;i<2;++i)
    {
        MeshLib::Location const location(_mesh.getID(),MeshLib::MeshItemType::Node,element.getNode(i)->getID());
        result.push_back(_local_to_global_index_map->getGlobalIndex(location,0,0));
        result.push_back(_local_to_global_index_map->getGlobalIndex(location,0,1));
        result.push_back(_local_to_global_index_map->getGlobalIndex(location,1,0));
    }
    return result;
}
void Beam2DProcess::assembleConcreteProcess(double,double,
    std::vector<GlobalVector*> const&,std::vector<GlobalVector*> const&,
    int,GlobalMatrix&,GlobalMatrix& K,GlobalVector&)
{
    if(data_.axial_capacity>0.) { OGS_FATAL("Plastic beams require Newton assembly."); }
    for(auto const* element:_mesh.getElements())
    {
        auto const ids=indices(*element);
        K.add({ids,ids},assemblers_[element->getID()].stiffness());
    }
}
void Beam2DProcess::assembleWithJacobianConcreteProcess(double,double,
    std::vector<GlobalVector*> const& x,std::vector<GlobalVector*> const&,
    int const process_id,GlobalVector& b,GlobalMatrix& jac)
{
    for(auto const* element:_mesh.getElements())
    {
        auto const ids=indices(*element);
        auto const& K=assemblers_[element->getID()].stiffness();
        Eigen::Vector<double,6> const q=MathLib::toVector(x[process_id]->get(ids));
        if(data_.axial_capacity>0.)
        {
            auto const response=plastic_assemblers_[element->getID()].trial(q,accepted_sections_[element->getID()]);
            Eigen::Vector<double,6> const rhs=-response.force;
            b.add(ids,rhs);jac.add({ids,ids},response.tangent);
        }
        else
        {
            Eigen::Vector<double,6> const rhs=-K*q;
            b.add(ids,rhs);jac.add({ids,ids},K);
        }
    }
}
void Beam2DProcess::computeSecondaryVariableConcrete(double,double,
    std::vector<GlobalVector*> const& x,GlobalVector const&,int const process_id)
{
    for(auto const* element:_mesh.getElements())
    {
        auto const ids=indices(*element);
        Eigen::Vector<double,6> const q=MathLib::toVector(x[process_id]->get(ids));
        auto force=assemblers_[element->getID()].localForces(q);
        if(data_.axial_capacity>0.)
        {
            auto const id=element->getID();
            auto const response=plastic_assemblers_[id].trial(q,accepted_sections_[id]);
            force=plastic_assemblers_[id].localForces(response.force);
            for(int ip=0;ip<2;++ip)
            {
                Eigen::Map<Eigen::Vector2d>(section_forces_->data()+4*id+2*ip)=response.section_forces[ip];
                Eigen::Map<Eigen::Vector2d>(plastic_strains_->data()+4*id+2*ip)=response.states[ip].plastic_deformation;
            }
        }
        Eigen::Map<Eigen::Vector<double,6>>(forces_->data()+6*element->getID())=force;
    }
}
void Beam2DProcess::postAcceptedTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,double,int const process_id)
{
    if(data_.axial_capacity==0.) { return; }
    for(auto const* element:_mesh.getElements())
    {
        auto const id=element->getID();auto const ids=indices(*element);
        Eigen::Vector<double,6> const q=MathLib::toVector(x[process_id]->get(ids));
        accepted_sections_[id]=plastic_assemblers_[id].trial(q,accepted_sections_[id]).states;
    }
}
}
