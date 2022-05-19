//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TimeDerivativePostprocessor.h"

registerMooseObject("MooseApp", TimeDerivativePostprocessor);

InputParameters
TimeDerivativePostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription("Compute the time derivative of a postprocessor");
  params.addRequiredParam<PostprocessorName>("postprocessor", "The name of the postprocessor");
  return params;
}

TimeDerivativePostprocessor::TimeDerivativePostprocessor(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _pps_value(getPostprocessorValue("postprocessor")),
    _pps_value_old(getPostprocessorValueOld("postprocessor"))
{
}

void
TimeDerivativePostprocessor::initialize()
{
}

void
TimeDerivativePostprocessor::execute()
{
  _rate = (_pps_value - _pps_value_old) / _dt;
}

Real
TimeDerivativePostprocessor::getValue()
{
  return _rate;
}
