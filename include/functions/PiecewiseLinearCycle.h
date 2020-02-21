//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CycleTime.h"
#include "FunctionInterface.h"
#include "PiecewiseLinear.h"

/**
 * Function which provides a PiecewiseLinear continuous linear interpolation
 * of a provided (x,y) point data set.
 */
class PiecewiseLinearCycle : public PiecewiseLinear,
                             protected FunctionInterface {
public:
  static InputParameters validParams();
  PiecewiseLinearCycle(const InputParameters &parameters);

  /**
   * Get the value of the function (based on time only)
   * \param t The time
   * \param pt The point in space (x,y,z) (unused)
   * \return The value of the function at the specified time
   */
  Real value(Real t, const Point &pt) const override;

  /**
   * Get the time derivative of the function (based on time only)
   * \param t The time
   * \param pt The point in space (x,y,z) (unused)
   * \return The time derivative of the function at the specified time
   */
  Real timeDerivative(Real t, const Point &pt) const override;

  Real integral() const override;

  Real average() const override;

protected:
  const Function &_cycle_time_func;
};
