//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CumulativeValuePostprocessor.h"

/**
 * Computes the time integral of a postprocessor representing a rate.
 */
class TimeIntegralPostprocessor : public CumulativeValuePostprocessor {
public:
  static InputParameters validParams();

  TimeIntegralPostprocessor(const InputParameters &parameters);

  virtual void execute() override;
};
