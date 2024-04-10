//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADStokesStressDivergence.h"

registerMooseObject("DeerApp", ADStokesStressDivergence);

InputParameters
ADStokesStressDivergence::validParams()
{
  InputParameters params = ADVectorKernel::validParams();
  params.addClassDescription("The unstabilized stress diveregence kernel for Stokes flow");

  return params;
}

ADStokesStressDivergence::ADStokesStressDivergence(const InputParameters & parameters)
  : ADVectorKernel(parameters),
    _deviatoric_stress(getADMaterialPropertyByName<RankTwoTensor>("deviatoric_stress"))
{
}

ADReal
ADStokesStressDivergence::computeQpResidual()
{
  return _deviatoric_stress[_qp].contract(_grad_test[_i][_qp]);
}
