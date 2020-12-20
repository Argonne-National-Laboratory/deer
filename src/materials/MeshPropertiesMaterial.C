//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshPropertiesMaterial.h"

registerMooseObject("DeerApp", MeshPropertiesMaterial);

InputParameters MeshPropertiesMaterial::validParams() {
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Compute the Mesh hmin, and hmax value as a material property");

  return params;
}

MeshPropertiesMaterial::MeshPropertiesMaterial(
    const InputParameters &parameters)
    : Material(parameters), _h_min(declareProperty<Real>("h_min")),
      _h_max(declareProperty<Real>("h_max")),
      _jxw(declareProperty<Real>("jxw")) {}

void MeshPropertiesMaterial::computeQpProperties() {
  _h_min[_qp] = _current_elem->hmin();
  _h_max[_qp] = _current_elem->hmax();
  _jxw[_qp] = _JxW[_qp];
}
