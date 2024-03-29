//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PiecewiseLinearCycle.h"

registerMooseObject("DeerApp", PiecewiseLinearCycle);

InputParameters
PiecewiseLinearCycle::validParams()
{
  InputParameters params = PiecewiseLinear::validParams();
  params.addRequiredParam<FunctionName>("cycle_time_func", "The CycleTime function name");
  params.addClassDescription("linearly interpolate a function cyclicly");
  return params;
}

PiecewiseLinearCycle::PiecewiseLinearCycle(const InputParameters & parameters)
  : PiecewiseLinear(parameters),
    FunctionInterface(this),
    _cycle_time_func(getFunction("cycle_time_func"))
{
}

Real
PiecewiseLinearCycle::value(Real t, const Point & p) const
{
  // Real c_fraction = _cycle_time_func::value(t, p);
  Real func_value;
  if (_has_axis)
  {
    func_value = _linear_interp->sample(p(_axis));
  }
  else
  {
    func_value = _linear_interp->sample(_cycle_time_func.value(t, p));
  }
  return _scale_factor * func_value;
}

Real
PiecewiseLinearCycle::timeDerivative(Real t, const Point & p) const
{
  Real func_value;
  if (_has_axis)
  {
    func_value = _linear_interp->sampleDerivative(p(_axis));
  }
  else
  {
    func_value = _linear_interp->sampleDerivative(_cycle_time_func.value(t, p));
  }
  std::cout << "PiecewiseLinearCycle: func_value " << func_value << "\n";
  return _scale_factor * func_value;
}

Real
PiecewiseLinearCycle::integral() const
{
  return _scale_factor * _linear_interp->integrate();
}

Real
PiecewiseLinearCycle::average() const
{
  return integral() /
         (_linear_interp->domain(_linear_interp->getSampleSize() - 1) - _linear_interp->domain(0));
}
