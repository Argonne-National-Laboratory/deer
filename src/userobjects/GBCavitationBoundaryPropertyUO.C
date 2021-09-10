//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBCavitationBoundaryPropertyUO.h"
#include "InterfaceValueTools.h"

registerMooseObject("DeerApp", GBCavitationBoundaryPropertyUO);

InputParameters GBCavitationBoundaryPropertyUO::validParams() {
  InputParameters params = InterfaceUserObject::validParams();
  params.addRequiredParam<FileName>("file_name",
                                    "the filename containing the different "
                                    "material properties for each GB pair");
  params.addClassDescription("Interface user obejcet storing and returing "
                             "properties for each grain boundary");
  return params;
}

GBCavitationBoundaryPropertyUO::GBCavitationBoundaryPropertyUO(
    const InputParameters &parameters)
    : InterfaceUserObject(parameters),
      _file_name(getParam<FileName>("file_name")) {}

void GBCavitationBoundaryPropertyUO::initialSetup() {
  ReadGBPropertyFile();
  constructGBPairsMap();
}

void GBCavitationBoundaryPropertyUO::constructGBPairsMap() {
  /* define the boundary map*/
  // build the global element,side,boundary_id list
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
      elem_side_bid = _mesh.buildSideList();

  // retrieve on which boundary this UO operates on
  std::set<BoundaryID> boundaryList = boundaryIDs();

  // initialize the map_values looping over all the element and sides
  for (unsigned int i = 0; i < elem_side_bid.size(); i++) {
    // check if this element side paris is part of the boundary this UO operates
    if (boundaryList.find(std::get<2>(elem_side_bid[i])) !=
        boundaryList.end()) {

      // I would like to avoid this but apparently we need to query the mesh to
      // get the proper nieghhbor at his stage.
      const Elem *current_elem = _mesh.elemPtr(std::get<0>(elem_side_bid[i]));
      const Elem *neighbor_elem =
          current_elem->neighbor_ptr(std::get<1>(elem_side_bid[i]));
      // construct the elemnt side pair
      std::pair<dof_id_type, unsigned int> elem_side_pair = std::make_pair(
          std::get<0>(elem_side_bid[i]), std::get<1>(elem_side_bid[i]));

      // query element and neighbor to get their subdomain ids to create a pair
      std::pair<SubdomainID, SubdomainID> subdomain_pair = std::make_pair(
          current_elem->subdomain_id(), neighbor_elem->subdomain_id());

      // add entry to the value elmeside_gbpair_map, will reuse this map in the
      // get mp values map
      _elmeside_gbpair_map[elem_side_pair] = subdomain_pair;
    }
  }
}

void GBCavitationBoundaryPropertyUO::ReadGBPropertyFile() {

  // open file file
  std::ifstream inFile(_file_name.c_str());
  if (!inFile)
    mooseError("Can't open ", _file_name);

  // Skip first 4 lines
  for (unsigned int i = 0; i < 4; ++i)
    inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  SubdomainID gid1, gid2;
  Real FN_NI, Nmax_NI, thickness, a_0, b_0, sigma_0, D_GB, beta_exponent,
      n_exponent, E_GB, G_GB, eta_sliding, psi_degree;

  while (inFile >> gid1 >> gid2 >> FN_NI >> Nmax_NI >> thickness >> a_0 >>
         b_0 >> sigma_0 >> D_GB >> beta_exponent >> n_exponent >> E_GB >>
         G_GB >> eta_sliding >> psi_degree) {

    std::pair<SubdomainID, SubdomainID> gb_pair = std::make_pair(gid1, gid2);
    std::map<std::string, Real> boundary_property;

    boundary_property["FN_NI"] = FN_NI;
    boundary_property["Nmax_NI"] = Nmax_NI;
    boundary_property["thickness"] = thickness;
    boundary_property["a0"] = a_0;
    boundary_property["b0"] = b_0;
    boundary_property["sigma_0"] = sigma_0;
    boundary_property["D_GB"] = D_GB;
    boundary_property["beta_exponent"] = beta_exponent;
    boundary_property["n_exponent"] = n_exponent;
    boundary_property["E_GB"] = E_GB;
    boundary_property["G_GB"] = G_GB;
    boundary_property["eta_sliding"] = eta_sliding;
    boundary_property["psi_degree"] = psi_degree;

    gb_property_map[gb_pair] = boundary_property;
  }
}

std::map<std::string, Real>
GBCavitationBoundaryPropertyUO::getPropertyMap(const dof_id_type elem_id,
                                               const unsigned int side) const {

  // find the subdomain pairs give an lement and a side
  std::pair<SubdomainID, SubdomainID> gb_pairs;
  auto gbpair_it = _elmeside_gbpair_map.find(std::make_pair(elem_id, side));
  if (gbpair_it != _elmeside_gbpair_map.end())
    gb_pairs = gbpair_it->second;
  else
    mooseError("GBCavitationBoundaryPropertyUO: can't find the given "
               "elem_id and side: " +
               std::to_string(elem_id) + ", " + std::to_string(side));

  // retrive properties give a subdomain pairs
  auto property_map = gb_property_map.find(gb_pairs);
  if (property_map != gb_property_map.end())
    return property_map->second;
  else
    mooseError("GBCavitationBoundaryPropertyUO: can't find the given "
               "boundary between block " +
               std::to_string(std::get<0>(gb_pairs)) + " and " +
               std::to_string(std::get<1>(gb_pairs)));
}
