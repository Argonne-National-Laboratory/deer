//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ScalarKernel.h"

class HomogenizationConstraintScalarKernel : public ScalarKernel
{
 public:
  static InputParameters validParams();

  HomogenizationConstraintScalarKernel(const InputParameters & parameters);

  virtual void reinit();
  virtual void computeResidual();
  virtual void computeJacobian();

 protected:
  unsigned int _h;

};
