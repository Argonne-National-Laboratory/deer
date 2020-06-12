//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NullScalarKernel.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableScalar.h"
#include "Function.h"

registerMooseObject("DeerApp", NullScalarKernel);

InputParameters
NullScalarKernel::validParams()
{
  InputParameters params = ScalarKernel::validParams();

  return params;
}

NullScalarKernel::NullScalarKernel(const InputParameters & parameters)
  : ScalarKernel(parameters)
{
}

void
NullScalarKernel::reinit()
{
}

void
NullScalarKernel::computeResidual()
{
  return;
}

void
NullScalarKernel::computeJacobian()
{
  // Amusingly you need an explicit zero or PETSC whines
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  for (_i = 0; _i < ke.m(); _i++)
    ke(_i, _i) += 0.0;
}

