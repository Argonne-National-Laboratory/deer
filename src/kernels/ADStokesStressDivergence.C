//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADStokesStressDivergence.h"

registerMooseObject("MooseApp", ADStokesStressDivergence);

InputParameters
ADStokesStressDivergence::validParams()
{
  InputParameters params = ADVectorKernel::validParams();
  params.addClassDescription("The unstabilized stress diveregence kernel for Stokes flow");

  return params;
}

ADStokesStressDivergence::ADStokesStressDivergence(const InputParameters & parameters)
  : ADVectorKernel(parameters), _stress(getADMaterialPropertyByName<RankTwoTensor>("stress"))
{
}

ADReal
ADStokesStressDivergence::computeQpResidual()
{
  return _stress[_qp].contract(_grad_phi[_qp]);
}