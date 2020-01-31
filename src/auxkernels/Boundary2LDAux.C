//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "Boundary2LDAux.h"

registerMooseObject("DeerApp", Boundary2LDAux);

defineLegacyParams(Boundary2LDAux);

InputParameters Boundary2LDAux::validParams() {
  InputParameters params = AuxKernel::validParams();
  params.addParam<UserObjectName>(
      "map2LDelem_uo_name", "the name of the user object containing the "
                            "map between boundary and LD elements and data");
  params.addClassDescription(
      "An aux kernel saving boundary data on   lower dimensional elements");

  return params;
}

Boundary2LDAux::Boundary2LDAux(const InputParameters &parameters)
    : AuxKernel(parameters),
      _map_2_ld_elem_uo(getUserObject<Map2LDInterfaceMaterialPropertyUO>(
          "map2LDelem_uo_name")) {}

Real Boundary2LDAux::computeValue() {

  return _map_2_ld_elem_uo.getValueForLD(_current_elem->id(), _qp);
}
