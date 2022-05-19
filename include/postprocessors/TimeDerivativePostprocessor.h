//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

class TimeDerivativePostprocessor;

/**
 * Computes the time derivative of Postprocessor, i.e. (v-v_old)/dt
 */

class TimeDerivativePostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  TimeDerivativePostprocessor(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;

protected:
  /// cumulative sum of the post-processor value
  Real _rate;

  /// current post-processor value
  const PostprocessorValue & _pps_value;

  /// old post-processor value
  const PostprocessorValue & _pps_value_old;
};
