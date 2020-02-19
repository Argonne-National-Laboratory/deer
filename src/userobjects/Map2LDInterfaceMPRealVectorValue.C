//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Map2LDInterfaceMPRealVectorValue.h"
#include "MooseMesh.h"
registerMooseObject("DeerApp", Map2LDInterfaceMPRealVectorValue);

InputParameters Map2LDInterfaceMPRealVectorValue::validParams() {
  InputParameters params = Map2LDInterfaceMaterialPropertyUO::validParams();
  params.addRequiredParam<unsigned int>("component",
                                        "the copmponent of the real vector");
  params.addClassDescription(
      "The userobject to output a componet of a RealVectorValue "
      "material property on LD elements");
  return params;
}

Map2LDInterfaceMPRealVectorValue::Map2LDInterfaceMPRealVectorValue(
    const InputParameters &parameters)
    : Map2LDInterfaceMaterialPropertyUO(parameters),
      _mp_value(getMaterialPropertyByName<RealVectorValue>(
          getParam<MaterialPropertyName>("mp_name"))),
      _component(getParam<unsigned int>("component")) {}

Map2LDInterfaceMPRealVectorValue::~Map2LDInterfaceMPRealVectorValue() {}
