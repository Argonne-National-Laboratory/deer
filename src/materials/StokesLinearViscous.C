//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "StokesLinearViscous.h"

registerMooseObject("DeerApp", StokesLinearViscous);

InputParameters
StokesLinearViscous::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription("Calculate the stress for Stokes flow using the simple linear model.");

  params.addRequiredParam<MaterialPropertyName>("mu", "Viscosity");

  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

StokesLinearViscous::StokesLinearViscous(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _stress(declareADProperty<RankTwoTensor>("stress")),
    _strain_rate(getADMaterialPropertyByName<RankTwoTensor>("strain_rate")),
    _mu(getADMaterialProperty<Real>("mu"))
{
}

void
StokesLinearViscous::computeQpProperties()
{
  _stress[_qp] = _strain_rate[_qp] * _mu[_qp];
}