// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <fstream>
#include <memory>
#include "MeshLib/IO/readMeshFromFile.h"

namespace ProcessLib::HydroMechanics
{
/// Preprocessing only: derive upper exterior edges from the unchanged bulk mesh.
/// Head-dependent wet/dry classification remains in the runtime boundary.
inline void prepareAutomaticReservoirSurfaces(
    boost::property_tree::ptree& root,
    boost::property_tree::ptree& process,
    std::filesystem::path const& project_directory,
    std::filesystem::path const& output_directory)
{
    using namespace InitialGroundwaterDetail;
    auto controls=process.get_child_optional("hydraulic_controls");
    if (!controls) { return; }
    bool requested=false;
    for (auto const& [tag,c]:*controls)
    {
        if (c.count("automatic_pressure_mesh")) { fail("automatic_pressure_mesh is an internal reserved entry"); }
        requested|=c.count("exterior_surface")!=0;
    }
    if (!requested) { return; }
    if (process.get<std::string>("coupling_scheme","monolithic")!="monolithic" ||
        process.get<std::string>("type")!="HYDRO_MECHANICS")
    { fail("automatic reservoir surfaces require monolithic HM"); }
    auto const bulk_name=process.get<std::string>("mesh",
        std::filesystem::path(root.get_child("meshes").front().second.data()).stem().string());
    std::filesystem::path bulk_path;
    std::map<std::string,std::filesystem::path> paths;
    for (auto const& [tag,m]:root.get_child("meshes"))
    {
        if (tag!="mesh") { continue; }
        auto const path=project_directory/m.data();
        paths.emplace(path.stem().string(),path);
        if (path.stem()==bulk_name) { bulk_path=path; }
    }
    if (bulk_path.empty()) { fail("automatic reservoir bulk mesh is not declared"); }
    std::unique_ptr<MeshLib::Mesh> bulk(MeshLib::IO::readMeshFromFile(bulk_path.string()));
    if (!bulk || bulk->getDimension()!=2) { fail("automatic reservoir surfaces currently require a 2D bulk mesh"); }
    struct Edge { std::size_t element; std::vector<std::size_t> nodes; };
    std::map<std::pair<std::size_t,std::size_t>,std::vector<Edge>> edges;
    for (auto const* e:bulk->getElements())
    {
        if (e->getDimension()!=2) { fail("automatic reservoir bulk mesh must contain only 2D elements"); }
        for (unsigned i=0;i<e->getNumberOfEdges();++i)
        {
            std::unique_ptr<MeshLib::Element const> edge(e->getEdge(i));
            if (edge->getNumberOfNodes()!=2 && edge->getNumberOfNodes()!=3)
            { fail("automatic reservoir edges support linear/quadratic geometry only"); }
            Edge item{e->getID(),{}};
            for (auto const* n:edge->nodes()) { item.nodes.push_back(n->getID()); }
            auto const pair=std::minmax(item.nodes[0],item.nodes[1]);
            edges[{pair.first,pair.second}].push_back(std::move(item));
        }
    }
    // Explicit far-field controls own shared pressure nodes (e.g. top/right corner).
    std::set<std::size_t> excluded;
    for (auto const& [tag,c]:*controls)
    {
        if (c.count("exterior_surface")) { continue; }
        auto const name=c.get<std::string>("mesh");
        if (!paths.contains(name)) { fail("unknown explicit hydraulic boundary mesh " + name); }
        std::unique_ptr<MeshLib::Mesh> boundary(MeshLib::IO::readMeshFromFile(paths.at(name).string()));
        if (!boundary) { fail("cannot read explicit hydraulic boundary mesh"); }
        auto const* ids=boundary->getProperties().getPropertyVector<std::size_t>("bulk_node_ids",MeshLib::MeshItemType::Node,1);
        if (!ids) { fail("explicit hydraulic boundaries need bulk_node_ids for automatic corner ownership"); }
        excluded.insert(ids->begin(),ids->end());
    }
    auto const out=std::filesystem::absolute(output_directory);
    std::filesystem::create_directories(out);
    auto write=[&](std::string const& name,std::vector<Edge> const& selected,bool nodal)
    {
        std::set<std::size_t> ids;
        for (auto const& e:selected) { ids.insert(e.nodes.begin(),e.nodes.end()); }
        if (nodal) { for (auto id:excluded) { ids.erase(id); } }
        if (ids.empty()) { fail("automatic reservoir selection has no applicable nodes"); }
        std::map<std::size_t,std::size_t> local;
        for (auto id:ids) { local[id]=local.size(); }
        auto const path=out/(name+".vtu");
        std::ofstream f(path); f<<std::setprecision(17);
        auto const cells=nodal?ids.size():selected.size();
        f<<"<?xml version=\"1.0\"?><VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"><UnstructuredGrid><Piece NumberOfPoints=\""<<ids.size()<<"\" NumberOfCells=\""<<cells<<"\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
        for (auto id:ids) { auto const& n=*bulk->getNode(id);f<<n[0]<<' '<<n[1]<<' '<<n[2]<<' '; }
        f<<"</DataArray></Points><PointData><DataArray type=\"UInt64\" Name=\"bulk_node_ids\" format=\"ascii\">";
        for (auto id:ids) { f<<id<<' '; }
        f<<"</DataArray></PointData><CellData>";
        if (!nodal)
        {
            f<<"<DataArray type=\"UInt64\" Name=\"bulk_element_ids\" format=\"ascii\">";
            for (auto const& e:selected) { f<<e.element<<' '; }
            f<<"</DataArray>";
        }
        f<<"</CellData><Cells><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">";
        if (nodal) { for (auto const& [id,i]:local) { f<<i<<' '; } }
        else { for (auto const& e:selected) { for (auto id:e.nodes) { f<<local.at(id)<<' '; } } }
        f<<"</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">";
        std::size_t offset=0;
        if (nodal) { for (std::size_t i=0;i<cells;++i) { f<<++offset<<' '; } }
        else { for (auto const& e:selected) { offset+=e.nodes.size();f<<offset<<' '; } }
        f<<"</DataArray><DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
        if (nodal) { for (std::size_t i=0;i<cells;++i) { f<<1<<' '; } }
        else { for (auto const& e:selected) { f<<(e.nodes.size()==3?21:3)<<' '; } }
        f<<"</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>";
        f.close();if (!f) { fail("cannot write automatic reservoir boundary"); }
        ptree entry;entry.data()=path.string();root.get_child("meshes").push_back({"mesh",entry});
        return name;
    };
    struct Assignment { std::set<std::size_t> nodes; std::vector<double> heads; };
    std::vector<Assignment> assignments;
    for (auto& [tag,c]:*controls)
    {
        auto selector=c.get_child_optional("exterior_surface");
        if (!selector) { continue; }
        if (c.count("mesh") || c.get<std::string>("type")!="ReservoirSeepage")
        { fail("exterior_surface replaces mesh and requires ReservoirSeepage"); }
        bool const supplied=selector->get<std::string>("type")=="Mesh";
        keys(*selector,supplied?std::initializer_list<std::string>{"type","mesh"}:std::initializer_list<std::string>{"type","region","divide_at_x"});
        for (auto const& [key,value]:*selector) { keys(value,{}); }
        auto const region=supplied?std::string("All"):selector->get<std::string>("region");
        if ((!supplied && selector->get<std::string>("type")!="Automatic") ||
            (region!="All" && region!="Left" && region!="Right") ||
            ((region=="All")==static_cast<bool>(selector->get_optional<double>("divide_at_x"))))
        { fail("automatic exterior region requires All without split, or Left/Right with divide_at_x"); }
        double const split=region=="All"?0.:finiteNumber(*selector,"divide_at_x");
        auto const id=c.get<std::string>("<xmlattr>.id");
        if (id.empty() || id.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")!=std::string::npos)
        { fail("automatic reservoir requires a simple boundary ID"); }
        auto const name="__hm_exterior_"+id;
        if (paths.contains(name) || paths.contains(name+"_pressure")) { fail("reserved automatic boundary mesh name collision"); }
        std::set<std::pair<std::size_t,std::size_t>> supplied_edges;
        std::string supplied_name;
        if (supplied)
        {
            supplied_name=selector->get<std::string>("mesh");
            if (!paths.contains(supplied_name)) { fail("unknown reservoir surface mesh " + supplied_name); }
            std::unique_ptr<MeshLib::Mesh> surface(MeshLib::IO::readMeshFromFile(paths.at(supplied_name).string()));
            if (!surface || surface->getDimension()!=1) { fail("reservoir surface mesh must contain edges"); }
            auto const* map=surface->getProperties().getPropertyVector<std::size_t>("bulk_node_ids",MeshLib::MeshItemType::Node,1);
            if (!map || map->size()!=surface->getNumberOfNodes()) { fail("reservoir surface needs bulk_node_ids"); }
            for (auto const* e:surface->getElements())
            {
                auto const a=(*map)[e->getNode(0)->getID()];auto const b=(*map)[e->getNode(1)->getID()];
                auto const key=std::minmax(a,b);
                if (!supplied_edges.insert({key.first,key.second}).second) { fail("duplicate reservoir surface edge"); }
            }
        }
        std::vector<Edge> selected;
        Assignment assignment{{},{finiteNumber(c,"initial_head")}};
        for (auto const& [phase_tag,phase]:process.get_child("phases"))
        {
            if (auto heads=phase.get_child_optional("heads"))
            { for (auto const& [head_tag,head]:*heads) { if (head.get<std::string>("<xmlattr>.ref")==id) { assignment.heads.push_back(finiteNumber(head,"target_head")); } } }
        }
        for (auto const& [key,items]:edges)
        {
            if (items.size()>2) { fail("non-manifold bulk edge in automatic reservoir detection"); }
            if (items.size()!=1) { continue; }
            if (supplied && !supplied_edges.contains(key)) { continue; }
            auto item=items.front();
            auto const& a=*bulk->getNode(item.nodes[0]);auto const& b=*bulk->getNode(item.nodes[1]);
            Eigen::Vector2d tangent(b[0]-a[0],b[1]-a[1]);
            if (tangent.norm()==0.) { fail("zero-length automatic exterior edge"); }
            Eigen::Vector2d normal(tangent[1],-tangent[0]);
            Eigen::Vector2d centre=Eigen::Vector2d::Zero();
            auto const* element=bulk->getElement(item.element);
            for (unsigned j=0;j<element->getNumberOfBaseNodes();++j)
            { auto const& n=*element->getNode(j);centre+=Eigen::Vector2d(n[0],n[1]); }
            centre/=element->getNumberOfBaseNodes();
            Eigen::Vector2d const midpoint(.5*(a[0]+b[0]),.5*(a[1]+b[1]));
            if (normal.dot(midpoint-centre)<0.)
            { normal=-normal;std::swap(item.nodes[0],item.nodes[1]); }
            // Upper exterior only: omit vertical far fields and downward base.
            if (normal[1]<=1e-12*normal.norm()) { continue; }
            if ((region=="Left" && midpoint[0]>=split) || (region=="Right" && midpoint[0]<split)) { continue; }
            selected.push_back(item);assignment.nodes.insert(item.nodes.begin(),item.nodes.end());
        }
        if (selected.empty() || (supplied && selected.size()!=supplied_edges.size()))
        { fail("reservoir surface must select nonempty upper exterior edges of the bulk mesh"); }
        for (auto const& previous:assignments)
        {
            for (auto node:assignment.nodes)
            {
                if (!previous.nodes.contains(node) || excluded.contains(node)) { continue; }
                double const y=(*bulk->getNode(node))[1];
                if (*std::max_element(previous.heads.begin(),previous.heads.end())>=y ||
                    *std::max_element(assignment.heads.begin(),assignment.heads.end())>=y)
                { fail("automatic reservoir regions share potentially wetted nodes; separate regions above all water levels"); }
            }
        }
        assignments.push_back(std::move(assignment));
        c.put("mesh",supplied?supplied_name:write(name,selected,false));
        c.put("automatic_pressure_mesh",write(name+"_pressure",selected,true));
        c.erase("exterior_surface");
    }
}
} // namespace ProcessLib::HydroMechanics
