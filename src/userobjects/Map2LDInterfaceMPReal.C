//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Map2LDInterfaceMPReal.h"
#include "MooseMesh.h"
registerMooseObject("DeerApp", Map2LDInterfaceMPReal);

defineLegacyParams(Map2LDInterfaceMPReal);

InputParameters Map2LDInterfaceMPReal::validParams() {
  InputParameters params = Map2LDInterfaceMaterialPropertyUO::validParams();
  params.addClassDescription(
      "The userobject to output Real material property on LD elements");
  return params;
}

Map2LDInterfaceMPReal::Map2LDInterfaceMPReal(const InputParameters &parameters)
    : Map2LDInterfaceMaterialPropertyUO(parameters),
      _mp_value(getMaterialPropertyByName<Real>(
          getParam<MaterialPropertyName>("mp_name"))) {}

Map2LDInterfaceMPReal::~Map2LDInterfaceMPReal() {}
