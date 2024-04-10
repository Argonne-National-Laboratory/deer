//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADStokesPressure.h"

registerMooseObject("DeerApp", ADStokesPressure);

InputParameters
ADStokesPressure::validParams()
{
  InputParameters params = ADVectorKernel::validParams();
  params.addClassDescription("The pressure part of the velocity kernel for Stokes flow, do not integrate by parts");

  params.addRequiredCoupledVar("pressure", "The pressure");

  return params;
}

ADStokesPressure::ADStokesPressure(const InputParameters & parameters)
  : ADVectorKernel(parameters),
    _grad_pressure(adCoupledGradient("pressure"))
{
}

ADReal
ADStokesPressure::computeQpResidual()
{
  return _test[_i][_qp] * _grad_pressure[_qp];
}
