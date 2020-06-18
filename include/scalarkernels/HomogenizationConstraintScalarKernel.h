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

#include "HomogenizationConstraintIntegral.h"
#include "HomogenizationConstraintKernel.h"

class HomogenizationConstraintScalarKernel : public ScalarKernel
{
 public:
  static InputParameters validParams();

  HomogenizationConstraintScalarKernel(const InputParameters & parameters);

  virtual void reinit();
  virtual void computeResidual();
  virtual void computeJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);

 protected:
  const bool _ld;
  unsigned int _ndisp;

  unsigned int _num_hvars;
  std::vector<unsigned int> _homogenization_nums;

  unsigned int _h;
  const HomogenizationConstraintIntegral & _integrator;

  // Useful Voigt stuff
  const HomogenizationConstants::index_list _indices;
};
