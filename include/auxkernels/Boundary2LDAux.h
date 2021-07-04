//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "AuxKernel.h"
#include "Map2LDInterfaceMaterialPropertyUO.h"

/**
 * A base class for the various Material related AuxKernal objects
 */
class Boundary2LDAux : public AuxKernel {
public:
  static InputParameters validParams();

  /**
   * Class constructor
   * @param parameters The input parameters for this object
   */
  Boundary2LDAux(const InputParameters &parameters);

protected:
  virtual Real computeValue() override;

  const Map2LDInterfaceMaterialPropertyUO &_map_2_ld_elem_uo;
};
