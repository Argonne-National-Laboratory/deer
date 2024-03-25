//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADStokesIncompressibility.h"

registerMooseObject("MooseApp", ADStokesIncompressibility);

InputParameters
ADStokesIncompressibility::validParams()
{
  InputParameters params = ADVectorKernel::validParams();
  params.addClassDescription("The incompressibility equation for Stokes flow");

  params.addRequiredCoupledVar("displacement", "The vector-valued displacement");

  return params;
}

ADStokesIncompressibility::ADStokesIncompressibility(const InputParameters & parameters)
  : ADVectorKernel(parameters), _grad_disp(adCoupledVectorGradient("displacement"))
{
}

ADReal
ADStokesIncompressibility::computeQpResidual()
{
  return (_grad_disp[_qp](0, 0) + _grad_disp[_qp](1, 1) + _grad_disp[_qp](2, 2)) * _phi[_qp];
}