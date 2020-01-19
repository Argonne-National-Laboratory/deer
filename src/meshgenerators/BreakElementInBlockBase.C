//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BreakElementInBlockBase.h"
#include "InputParameters.h"

defineLegacyParams(BreakElementInBlockBase);

InputParameters BreakElementInBlockBase::validParams() {
  InputParameters params = MeshGenerator::validParams();

  params.addParam<std::string>("interface_name", "interface",
                               "the name of the new interface. Cannot be used "
                               "whit `split_interface=true`");
  params.addRequiredParam<std::vector<subdomain_id_type>>(
      "blocks_id", "the blocks in which element should be broken");
  params.addClassDescription(
      "This is the base class used to split all elemnts in a block");
  params.addParam<bool>("split_interface", false,
                        "If true, it create a "
                        "different interface for each block pair.");
  return params;
}

BreakElementInBlockBase::BreakElementInBlockBase(
    const InputParameters &parameters)
    : MeshGenerator(parameters),
      _interface_name(getParam<std::string>("interface_name")),
      _split_interface(getParam<bool>("split_interface")) {
  std::vector<subdomain_id_type> blocks =
      getParam<std::vector<subdomain_id_type>>("blocks_id");

  for (unsigned int i = 0; i < blocks.size(); i++)
    _blocks_id.insert(blocks[i]);
}

std::string BreakElementInBlockBase::generateBoundaryName(
    MeshBase &mesh, const subdomain_id_type &masterBlockID,
    const subdomain_id_type &slaveBlockID) {
  std::string master_block_name = mesh.subdomain_name(masterBlockID);
  std::string slave_block_name = mesh.subdomain_name(slaveBlockID);
  if (master_block_name.empty())
    master_block_name = "Block" + std::to_string(masterBlockID);
  if (slave_block_name.empty())
    slave_block_name = "Block" + std::to_string(slaveBlockID);

  return master_block_name + "_" + slave_block_name;
}

boundary_id_type BreakElementInBlockBase::findFreeBoundaryId(MeshBase &mesh) {
  const std::set<boundary_id_type> &currentBoundaryIds =
      mesh.get_boundary_info().get_boundary_ids();
  bool freeBoundaryNotFound = true;
  boundary_id_type freeId;
  for (freeId = 0; freeId < std::numeric_limits<boundary_id_type>::max();
       freeId++) {
    if (currentBoundaryIds.count(freeId) == 0) {
      // bid is not in the set, boundaryID is free
      freeBoundaryNotFound = false;
      break;
    }
  }

  if (freeBoundaryNotFound)
    mooseError("Too many boundaries. Maximum limit exceeded!");

  return freeId;
}

void BreakElementInBlockBase::mapBoundaryIdAndBoundaryName(
    boundary_id_type &boundaryID, const std::string &boundaryName) {
  _bName_bID_set.insert(std::pair<std::string, int>(boundaryName, boundaryID));
}
