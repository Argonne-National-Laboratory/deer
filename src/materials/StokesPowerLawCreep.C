//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "StokesPowerLawCreep.h"

registerMooseObject("DeerApp", StokesPowerLawCreep);

InputParameters
StokesPowerLawCreep::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription("Calculate the stress for Stokes flow using a power law creep model");

  params.addRequiredParam<MaterialPropertyName>("A", "Creep prefactor");
  params.addRequiredParam<MaterialPropertyName>("n", "Creep exponent");

  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

StokesPowerLawCreep::StokesPowerLawCreep(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _stress(declareADProperty<RankTwoTensor>("deviatoric_stress")),
    _strain_rate(getADMaterialPropertyByName<RankTwoTensor>("strain_rate")),
    _A(getADMaterialProperty<Real>("A")),
    _n(getADMaterialProperty<Real>("n"))
{
}

void
StokesPowerLawCreep::computeQpProperties()
{
  auto A = _A[_qp];
  auto n = _n[_qp];
  auto dev_strain_rate = _strain_rate[_qp].deviatoric();

  auto effective_strain_rate = std::sqrt(2.0 / 3.0 * dev_strain_rate.contract(dev_strain_rate));
  auto effective_stress = std::pow(effective_strain_rate / A, 1.0 / n);
  
  _stress[_qp] = 2.0 / 3.0 * dev_strain_rate / (A * std::pow(effective_stress, n - 1.0));
}
