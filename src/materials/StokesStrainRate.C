//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "StokesStrainRate.h"

registerMooseObject("DeerApp", StokesStrainRate);

InputParameters
StokesStrainRate::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription("Calculate the strain rate as the symmetric gradient of the velocity");
  params.addRequiredCoupledVar("velocity", "The vector-valued velocity");
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

StokesStrainRate::StokesStrainRate(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _strain_rate(declareADProperty<RankTwoTensor>("strain_rate")),
    _grad_vel(adCoupledVectorGradient("velocity"))

{
}

void
StokesStrainRate::computeQpProperties()
{
  _strain_rate[_qp] = 0.5 * (_grad_vel[_qp] + _grad_vel[_qp].transpose());
}