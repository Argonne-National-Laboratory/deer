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

/**
 *
 */
class Map2LDInterfaceMPRealVectorValue
    : public Map2LDInterfaceMaterialPropertyUO {
public:
  static InputParameters validParams();

  Map2LDInterfaceMPRealVectorValue(const InputParameters &parameters);
  virtual ~Map2LDInterfaceMPRealVectorValue();

protected:
  virtual Real getMPValue(const unsigned int qp) const override {
    return _mp_value[qp](_component);
  };
  const MaterialProperty<RealVectorValue> &_mp_value;
  const unsigned int &_component;
};
