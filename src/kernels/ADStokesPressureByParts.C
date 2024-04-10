//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADStokesPressureByParts.h"

registerMooseObject("DeerApp", ADStokesPressureByParts);

InputParameters
ADStokesPressureByParts::validParams()
{
  InputParameters params = ADVectorKernel::validParams();
  params.addClassDescription("The pressure part of the velocity kernel for Stokes flow, do not integrate by parts");

  params.addRequiredCoupledVar("pressure", "The pressure");

  return params;
}

ADStokesPressureByParts::ADStokesPressureByParts(const InputParameters & parameters)
  : ADVectorKernel(parameters),
    _pressure(adCoupledValue("pressure"))
{
}

ADReal
ADStokesPressureByParts::computeQpResidual()
{
  return -_grad_test[_i][_qp].tr() * _pressure[_qp];
}
