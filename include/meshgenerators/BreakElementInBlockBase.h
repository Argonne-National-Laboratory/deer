//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"

// Forward declarations
class BreakElementInBlockBase;

template <> InputParameters validParams<BreakElementInBlockBase>();

/**
 *
 */
class BreakElementInBlockBase : public MeshGenerator {
public:
  static InputParameters validParams();

  BreakElementInBlockBase(const InputParameters &parameters);

protected:
  /// the file_name from whence this mesh came
  std::string _file_name;
  /// the name of the new interface
  std::string _interface_name;
  /// the flag to split the interface by block
  bool _split_interface;
  /// a vector of block ids on which this mesh generator operates on
  std::set<subdomain_id_type> _blocks_id;

  std::set<std::pair<std::string, BoundaryID>> _bName_bID_set;

  /// this method finds the first free boundary id
  BoundaryID findFreeBoundaryId(MeshBase &mesh);

  /// this method generate the boundary name by assembling subdomain names
  std::string generateBoundaryName(MeshBase &mesh,
                                   const subdomain_id_type & /*masterBlockID*/,
                                   const subdomain_id_type & /*slaveBlockID*/);

private:
  /// this method save the boundary name/id pair
  void mapBoundaryIdAndBoundaryName(boundary_id_type & /*boundaryID*/,
                                    const std::string & /*boundaryName*/);
};
