//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADVectorKernel.h"

class ADStokesPressure : public ADVectorKernel
{
public:
  static InputParameters validParams();

  ADStokesPressure(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADVariableGradient & _grad_pressure;
};
