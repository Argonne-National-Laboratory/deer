//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "BreakElementInBlockBase.h"

class BreakElementInBlock;

template <> InputParameters validParams<BreakElementInBlock>();

class BreakElementInBlock : public BreakElementInBlockBase {
public:
  static InputParameters validParams();

  BreakElementInBlock(const InputParameters &parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  std::unique_ptr<MeshBase> &_input;

private:
  /// generate the new boundary interface
  void addInterfaceBoundary(MeshBase &mesh);

  std::set<std::pair<subdomain_id_type, subdomain_id_type>>
      _neighboring_block_list;
  std::map<std::pair<subdomain_id_type, subdomain_id_type>,
           std::set<std::pair<dof_id_type, unsigned int>>>
      _new_boundary_sides_map;
};
