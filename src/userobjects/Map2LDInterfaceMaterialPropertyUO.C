//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Map2LDInterfaceMaterialPropertyUO.h"
#include "MooseMesh.h"

defineLegacyParams(Map2LDInterfaceMaterialPropertyUO);

InputParameters Map2LDInterfaceMaterialPropertyUO::validParams() {
  InputParameters params = Map2LDelem::validParams();
  params.addRequiredParam<MaterialPropertyName>("mp_name",
                                                "The material property name");
  params.addClassDescription(
      "The Base userobject to output material properties on LD elements. The "
      "output for all material types shold subclass this");
  return params;
}

Map2LDInterfaceMaterialPropertyUO::Map2LDInterfaceMaterialPropertyUO(
    const InputParameters &parameters)
    : Map2LDelem(parameters) {}

Map2LDInterfaceMaterialPropertyUO::~Map2LDInterfaceMaterialPropertyUO() {}

void Map2LDInterfaceMaterialPropertyUO::initialSetup() {
  Map2LDelem::initialSetup();
  // define the boundary map and retrieve element side and boundary_ID
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
      elem_side_bid = _mesh.buildSideList();

  // retrieve on which boundary this UO operates
  std::set<BoundaryID> boundaryList = boundaryIDs();

  // clear map values
  _map_values.clear();

  // initialize the map_values looping over all the element and sides
  for (unsigned int i = 0; i < elem_side_bid.size(); i++) {
    // check if this element side is part of the boundary, if so add element
    // side to the interface map
    if (boundaryList.find(std::get<2>(elem_side_bid[i])) !=
        boundaryList.end()) {
      // make pair
      std::pair<dof_id_type, unsigned int> elem_side_pair = std::make_pair(
          std::get<0>(elem_side_bid[i]), std::get<1>(elem_side_bid[i]));
      // initialize map elemenet
      std::vector<Real> var_values(0, 0);

      // add entry to the value map
      _map_values[elem_side_pair] = var_values;
    }
  }
}

void Map2LDInterfaceMaterialPropertyUO::execute() {
  // find the entry on the map
  auto it =
      _map_values.find(std::make_pair(_current_elem->id(), _current_side));
  if (it != _map_values.end()) {
    // insert two vector value for each qp
    auto &vec = _map_values[std::make_pair(_current_elem->id(), _current_side)];
    vec.resize(_qrule->n_points());

    // loop over qps and do stuff
    for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
      // compute average value at qp
      vec[qp] = getMPValue(qp);
  } else
    mooseError(
        "Map2LDInterfaceMaterialPropertyUO:: cannot find the required element "
        "and side");
}

Real Map2LDInterfaceMaterialPropertyUO::getQpValue(dof_id_type elem,
                                                   unsigned int side,
                                                   unsigned int qp) const {
  auto data = _map_values.find(std::make_pair(elem, side));
  if (data != _map_values.end())
    return data->second[qp];
  else
    mooseError("getMeanMatProp: can't find the given qp");
}

Real Map2LDInterfaceMaterialPropertyUO::getValueForLD(dof_id_type ld_elem,
                                                      unsigned int qp) const {

  std::pair<dof_id_type, unsigned int> elem_side = getLDNeighbor(ld_elem);

  return getQpValue(std::get<0>(elem_side), std::get<1>(elem_side), qp);
}
