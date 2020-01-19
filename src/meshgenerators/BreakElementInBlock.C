//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BreakElementInBlock.h"
#include "CastUniquePointer.h"

#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"

#include <typeinfo>

registerMooseObject("DeerApp", BreakElementInBlock);

defineLegacyParams(BreakElementInBlock);

InputParameters BreakElementInBlock::validParams() {
  InputParameters params = BreakElementInBlockBase::validParams();
  params.addRequiredParam<MeshGeneratorName>("input",
                                             "The mesh we want to modify");
  params.addClassDescription(
      "Break boundaries based on the subdomains to which their sides are "
      "attached. Naming convention for the new boundaries will be the old "
      "boundary name plus \"_to_\" plus the subdomain name. At the moment"
      "this only works on REPLICATED mesh");
  return params;
}

BreakElementInBlock::BreakElementInBlock(const InputParameters &parameters)
    : BreakElementInBlockBase(parameters), _input(getMesh("input")) {
  if (typeid(_input).name() == typeid(DistributedMesh).name())
    mooseError("BreakElementInBlock only works with ReplicatedMesh.");
}

std::unique_ptr<MeshBase> BreakElementInBlock::generate() {
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  const dof_id_type max_node_id = mesh->max_node_id();
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  const unique_id_type max_unique_id = mesh->parallel_max_unique_id();
#endif

  // initialize the node to element map
  std::map<dof_id_type, std::set<dof_id_type>> node_to_elem_map;
  for (const auto &elem : mesh->active_element_ptr_range()) {
    for (unsigned int n = 0; n < elem->n_nodes(); n++)
      node_to_elem_map[elem->node_id(n)].insert(elem->id());
  }

  // create the node_elem_to_replicate map
  std::map<dof_id_type, std::set<dof_id_type>> node_elem_to_replicate;
  for (auto node_it = node_to_elem_map.begin();
       node_it != node_to_elem_map.end(); ++node_it) {
    bool boundary_node = false;
    dof_id_type node_id = mesh->node_ptr(node_it->first)->id();
    for (auto elem_id = node_it->second.begin();
         elem_id != node_it->second.end(); elem_id++) {
      const Elem *current_elem = mesh->elem_ptr(*elem_id);
      if (_blocks_id.find(current_elem->subdomain_id()) == _blocks_id.end()) {
        boundary_node = true;
        break;
      }
    }
    if (boundary_node == false) {
      for (auto elem_id = node_it->second.begin();
           elem_id != node_it->second.end(); elem_id++)
        node_elem_to_replicate[node_id].insert(*elem_id);
      node_elem_to_replicate[node_id].erase(
          node_elem_to_replicate[node_id].begin());
    } else {
      for (auto elem_id = node_it->second.begin();
           elem_id != node_it->second.end(); elem_id++) {
        const Elem *current_elem = mesh->elem_ptr(*elem_id);
        if (_blocks_id.find(current_elem->subdomain_id()) != _blocks_id.end())
          node_elem_to_replicate[node_id].insert(*elem_id);
      }
    }
  }

  dof_id_type node_counter = max_node_id + 1;
  // loop over elements
  for (const auto &elem : mesh->active_element_ptr_range()) {
    // check if we need to replicate the nodes of this element
    if (_blocks_id.find(elem->subdomain_id()) != _blocks_id.end()) {
      dof_id_type elem_id = elem->id();
      for (unsigned int node_id = 0; node_id < elem->n_nodes(); ++node_id) {

        const Node *current_node = mesh->node_ptr(elem->node_id(node_id));

        // check if have to add this node in this element
        if (node_elem_to_replicate[current_node->id()].find(elem_id) !=
            node_elem_to_replicate[current_node->id()].end()) {

          // add new node
          Node *new_node = Node::build(*current_node, node_counter).release();

          new_node->processor_id() = elem->processor_id();
          mesh->add_node(new_node);
          // Add boundary info to the new node
          std::vector<boundary_id_type> node_boundary_ids;
          mesh->boundary_info->boundary_ids(current_node, node_boundary_ids);
          mesh->boundary_info->add_node(new_node, node_boundary_ids);

          elem->set_node(node_id) = new_node;
          node_counter++;
        }
      }
    }
  }

  // // delete unsed nodes
  // for (auto node_it = can_delete_node_map.begin();
  //      node_it != can_delete_node_map.end(); ++node_it)
  //   if (node_it->second == true)
  //     mesh->delete_node(mesh->node_ptr(node_it->first));

  // // compact nodes id
  // dof_id_type counter = 0;
  // for (auto node_ptr : mesh->node_ptr_range()) {
  //   if (node_ptr->valid_id() != true) {
  //     node_ptr->set_id() = counter;
  //     // node_ptr->set_unique_id() = counter;
  //   }
  //   counter++;
  // }
  // loop over all the elements ti generate the interface map
  for (const auto &elem : mesh->active_element_ptr_range()) {
    subdomain_id_type elem_subdomain = elem->subdomain_id();
    if (_blocks_id.find(elem_subdomain) != _blocks_id.end()) {
      dof_id_type elem_id = elem->id();
      for (auto neighbor : elem->neighbor_ptr_range()) {
        if (neighbor != nullptr) {
          // create subdomain pair
          subdomain_id_type neighbor_subdomain = neighbor->subdomain_id();
          std::pair<subdomain_id_type, subdomain_id_type> blocks_pair;
          if (elem_subdomain > neighbor_subdomain)
            blocks_pair = std::make_pair(neighbor_subdomain, elem_subdomain);
          else
            blocks_pair = std::make_pair(elem_subdomain, neighbor_subdomain);

          // find neighboring side
          std::pair<dof_id_type, unsigned int> elem_side_pair;
          dof_id_type neighbor_id = neighbor->id();
          if (elem_id < neighbor_id)
            elem_side_pair =
                std::make_pair(elem_id, elem->which_neighbor_am_i(neighbor));
          else
            elem_side_pair = std::make_pair(
                neighbor_id, neighbor->which_neighbor_am_i(elem));

          _new_boundary_sides_map[blocks_pair].insert(elem_side_pair);
        }
      }
    }
  }

  addInterfaceBoundary(*mesh);
  // mesh->prepare_for_use();
  return dynamic_pointer_cast<MeshBase>(mesh);
}

void BreakElementInBlock::addInterfaceBoundary(MeshBase &mesh) {
  BoundaryInfo &boundary_info = mesh.get_boundary_info();
  std::set<boundary_id_type> ids = boundary_info.get_boundary_ids();
  boundary_id_type boundaryID = *ids.rbegin() + 1;

  // Make sure the new is the same on every processor
  mesh.comm().max(boundaryID);

  if (_split_interface == false)
    boundary_info.sideset_name(boundaryID) = _interface_name;

  // loop over boundary sides
  for (auto &boundary_side_map : _new_boundary_sides_map) {

    if (_split_interface == true)
      boundary_info.sideset_name(boundaryID) = generateBoundaryName(
          mesh, boundary_side_map.first.first, boundary_side_map.first.second);

    // loop over all the side belonging to each pair and add it to the
    // proper interface
    for (auto &element_side_pair : boundary_side_map.second) {
      boundary_info.add_side(element_side_pair.first, element_side_pair.second,
                             boundaryID);
    }
    if (_split_interface == true)
      ++boundaryID;
  }
}
