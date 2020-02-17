//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CycleNumber.h"

#include "MathUtils.h"

registerMooseObject("DeerApp", CycleNumber);

InputParameters CycleNumber::validParams() {
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("cycle_period", "period of the cycle");
  params.addClassDescription("Function providing the cycle number given "
                             "simulation time and cycle period");
  return params;
}

CycleNumber::CycleNumber(const InputParameters &parameters)
    : Function(parameters), _cycle_period(getParam<Real>("cycle_period")) {}

Real CycleNumber::value(Real /*t*/, const Point & /*p*/) const {
  Real intpart;
  Real fractpart = std::modf(_t / _cycle_period, &intpart);
  return intpart;
}
