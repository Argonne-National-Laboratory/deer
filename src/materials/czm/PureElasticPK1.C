//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PureElasticPK1.h"

registerMooseObject("DeerApp", PureElasticPK1);

InputParameters PureElasticPK1::validParams() {
  InputParameters params = CZMMaterialBasePK1::validParams();
  params.addClassDescription("3D Coupled (3DC) cohesive law of Salehani and Irani with no damage");
  params.addParam<Real>("E", 1e3, "opening stiffness");
  params.addParam<Real>("G", 1e2, "shear stiffness");
  return params;
}

PureElasticPK1::PureElasticPK1(const InputParameters & parameters)
  : CZMMaterialBasePK1(parameters),
    _K(std::vector<Real>{getParam<Real>("E"), getParam<Real>("G"), getParam<Real>("G")})
{
}

void
PureElasticPK1::computeTractionIncrementAndDerivatives()
{
  _traction_inc[_qp] = _K * _displacement_jump_inc[_qp];
  _dtraction_djump[_qp] = _K;
}
