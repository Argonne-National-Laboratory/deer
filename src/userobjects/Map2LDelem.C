//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Map2LDelem.h"
#include "MooseMesh.h"
registerMooseObject("DeerApp", Map2LDelem);

defineLegacyParams(Map2LDelem);

InputParameters Map2LDelem::validParams() {
  InputParameters params = InterfaceUserObject::validParams();
  params.addRequiredParam<std::vector<SubdomainName>>(
      "ld_block_names", "The name of the lower dimensional element blocks");
  return params;
}

Map2LDelem::Map2LDelem(const InputParameters &parameters)
    : InterfaceUserObject(parameters),
      ld_block_ids(_mesh.getSubdomainIDs(
          getParam<std::vector<SubdomainName>>("ld_block_names"))) {}

Map2LDelem::~Map2LDelem() {}

void Map2LDelem::initialSetup() {
  // define the boundary map nad retrieve element side and boundary_ID
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
      elem_side_bid = _mesh.buildSideList();

  const std::map<dof_id_type, std::vector<dof_id_type>> &n_2_elem_map =
      _mesh.nodeToElemMap();

  // retrieve on which boudnary this UO operates
  std::set<BoundaryID> boundaryList = boundaryIDs();

  // clear map values
  _map_LD_with_elem_side.clear();

  // initialize the map_values looping over all the element and sides
  for (unsigned int i = 0; i < elem_side_bid.size(); i++) {
    dof_id_type LD_neighbor_elem_id = -1;
    // if this element side is part of the boundary then add elements to the map
    if (boundaryList.find(std::get<2>(elem_side_bid[i])) !=
        boundaryList.end()) {

      // need to find the LD element
      // start from the current element
      const Elem *elem_ptr = _mesh.elemPtr(std::get<0>(elem_side_bid[i]));
      // find all the nodes generatign the current side (local)
      std::vector<unsigned int> nodes_on_boundary_side =
          elem_ptr->nodes_on_side(std::get<1>(elem_side_bid[i]));
      // generate a set containing the global node ids of the current side
      std::unordered_set<unsigned int> nodes_on_current_side;
      for (unsigned int k = 0; k < nodes_on_boundary_side.size(); k++)
        nodes_on_current_side.insert(
            elem_ptr->node_id(nodes_on_boundary_side[k]));

      // find global node id of current local node
      unsigned int global_node_id =
          elem_ptr->node_id(nodes_on_boundary_side[0]);
      // list all the lements connected to this node
      auto list_of_connected_elem = n_2_elem_map.find(global_node_id);
      if (list_of_connected_elem !=
          n_2_elem_map.end()) // if there are elements (this might be removed)
      {
        bool is_LD_neighbor = true;
        // loop over all the connceted elements
        for (auto elem_id : list_of_connected_elem->second) {
          is_LD_neighbor = true;
          const Elem *neighbor_ptr = _mesh.elemPtr(elem_id);
          LD_neighbor_elem_id = neighbor_ptr->id();

          // if the current eleemnt is a LD element
          if (neighbor_ptr->dim() == (elem_ptr->dim() - 1)) {
            // check if this LD element has all the node in common with the face
            unsigned int n_elem_nodes = neighbor_ptr->n_nodes();

            // we loop over all its node to check that all of them belong to the
            // current side
            for (unsigned int k = 0; k < n_elem_nodes; k++)
              if (nodes_on_current_side.find(neighbor_ptr->node_id(k)) ==
                  nodes_on_current_side.end()) {
                is_LD_neighbor = false;
                break; // if we can't find a node we break the loop
              }
          } else
            is_LD_neighbor = false;

          if (is_LD_neighbor) // found a neighbor we need to add the element to
                              // the map
          {
            // we need to add the elelment on the final map

            // make pair
            std::pair<dof_id_type, unsigned int> elem_side_pair =
                std::make_pair(std::get<0>(elem_side_bid[i]),
                               std::get<1>(elem_side_bid[i]));
            // initialize map elemenet
            // add entry to the <element,side> LD neighbor map and to the
            // reverse map

            _map_LD_with_elem_side[LD_neighbor_elem_id] = elem_side_pair;
            break;
          }
        }                    // end of loop over connected element
        if (!is_LD_neighbor) // we couldn't find a neighbor
          mooseError("can't find LD neighbor for element " +
                     std::to_string(elem_ptr->id()) + " on side " +
                     std::to_string(std::get<1>(elem_side_bid[i])));
      } else
        mooseError("can't find any element connected to this node");
    }
  }
}

std::pair<dof_id_type, unsigned int>
Map2LDelem::getLDNeighbor(dof_id_type ld_elem) const {
  auto elem_side_pair = _map_LD_with_elem_side.find(ld_elem);
  if (elem_side_pair != _map_LD_with_elem_side.end()) {
    return elem_side_pair->second;
  } else
    mooseError("getLDNeighbor: can't find the given LD element");
}
