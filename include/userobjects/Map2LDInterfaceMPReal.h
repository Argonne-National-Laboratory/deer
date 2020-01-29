//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Map2LDInterfaceMaterialPropertyUO.h"

class Map2LDInterfaceMPReal;

template <> InputParameters validParams<Map2LDInterfaceMPReal>();

/**
 *
 */
class Map2LDInterfaceMPReal : public Map2LDInterfaceMaterialPropertyUO {
public:
  static InputParameters validParams();

  Map2LDInterfaceMPReal(const InputParameters &parameters);
  virtual ~Map2LDInterfaceMPReal();

protected:
  virtual Real getMPValue(const unsigned int qp) const override {
    return _mp_value[qp];
  };
  const MaterialProperty<Real> &_mp_value;
};
