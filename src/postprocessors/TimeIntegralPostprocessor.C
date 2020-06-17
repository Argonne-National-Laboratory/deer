//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TimeIntegralPostprocessor.h"

registerMooseObject("DeerApp", TimeIntegralPostprocessor);

InputParameters TimeIntegralPostprocessor::validParams() {
  InputParameters params = CumulativeValuePostprocessor::validParams();
  params.addClassDescription("Compute the time integral of a postprocessor "
                             "representing a rate.");
  return params;
}

TimeIntegralPostprocessor::TimeIntegralPostprocessor(
    const InputParameters &parameters)
    : CumulativeValuePostprocessor(parameters) {}

void TimeIntegralPostprocessor::execute() {
  _sum = _sum_old + _pps_value * _dt;
}
