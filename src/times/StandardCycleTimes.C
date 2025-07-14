//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StandardCycleTimes.h"

registerMooseObject("MooseApp", StandardCycleTimes);

InputParameters
StandardCycleTimes::validParams()
{
  InputParameters params = Times::validParams();
  params.addClassDescription("Times set per our standard creep-fatigue cycle");
  params.addRequiredParam<Real>("load", "Load time");
  params.addRequiredParam<size_t>("number_load", "Number of times in each loading segment");
  params.addRequiredParam<Real>("tension_hold", "Tension hold time");
  params.addRequiredParam<size_t>("number_tension_hold", "Number of times in each tension hold");
  params.addRequiredParam<Real>("compression_hold", "Compression hold time");
  params.addRequiredParam<size_t>("number_compression_hold", "Number of times in each compression hold");
  params.addRequiredParam<size_t>("cycles", "Number of load cycles total");

  // Times are known for all processes already
  params.set<bool>("auto_broadcast") = false;

  return params;
}

StandardCycleTimes::StandardCycleTimes(const InputParameters & parameters) : Times(parameters), 
    _load(getParam<Real>("load")), _tension_hold(getParam<Real>("tension_hold")), _compression_hold(getParam<Real>("compression_hold")), 
    _number_load(getParam<size_t>("number_load")), _number_tension_hold(getParam<size_t>("number_tension_hold")),
    _number_compression_hold(getParam<size_t>("number_compression_hold")),
    _cycles(getParam<size_t>("cycles"))
{
  std::vector<Real> times = {0.0};
  for (size_t i = 0; i < _cycles; i++)
  {
    for (size_t j = 0; j < _number_load; j++)
      times.push_back(times.back() + _load / _number_load);
    for (size_t j = 0; j < _number_tension_hold; j++)
      times.push_back(times.back() + _tension_hold / _number_tension_hold);
    for (size_t j = 0; j < _number_load; j++)
      times.push_back(times.back() + _load / _number_load);
    for (size_t j = 0; j < _number_load; j++)
      times.push_back(times.back() + _load / _number_load);
    for (size_t j = 0; j < _number_compression_hold; j++)
      times.push_back(times.back() + _compression_hold / _number_compression_hold);
    for (size_t j = 0; j < _number_load; j++)
      times.push_back(times.back() + _load / _number_load);
  }
  _times = times;
}
