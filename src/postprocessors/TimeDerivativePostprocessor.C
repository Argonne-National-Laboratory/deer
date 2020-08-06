//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TimeDerivativePostprocessor.h"

registerMooseObject("DeerApp", TimeDerivativePostprocessor);

defineLegacyParams(TimeDerivativePostprocessor);

InputParameters TimeDerivativePostprocessor::validParams() {
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<PostprocessorName>("value", "First value");
  params.addClassDescription("Computes the sum of two postprocessors");
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
  return params;
}

TimeDerivativePostprocessor::TimeDerivativePostprocessor(
    const InputParameters &parameters)
    : GeneralPostprocessor(parameters), _value(getPostprocessorValue("value")),
      _value_old(getPostprocessorValueOld("value")) {}

void TimeDerivativePostprocessor::initialize() {}

void TimeDerivativePostprocessor::execute() {}

PostprocessorValue TimeDerivativePostprocessor::getValue() {
  return (_value - _value_old) / _dt;
}
