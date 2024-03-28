//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADStokesIncompressibility.h"

registerMooseObject("DeerApp", ADStokesIncompressibility);

InputParameters
ADStokesIncompressibility::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The incompressibility equation for Stokes flow");

  params.addRequiredCoupledVar("velocity", "The vector-valued velocity");

  return params;
}

ADStokesIncompressibility::ADStokesIncompressibility(const InputParameters & parameters)
  : ADKernel(parameters), _grad_vel(adCoupledVectorGradient("velocity"))
{
}

ADReal
ADStokesIncompressibility::computeQpResidual()
{
  return -(_grad_vel[_qp](0, 0) + _grad_vel[_qp](1, 1) + _grad_vel[_qp](2, 2)) * _test[_i][_qp];
}