//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceAreaPostprocessor.h"

registerMooseObject("DeerApp", InterfaceAreaPostprocessor);

defineLegacyParams(InterfaceAreaPostprocessor);

InputParameters InterfaceAreaPostprocessor::validParams() {
  InputParameters params = InterfaceIntegralPostprocessor::validParams();

  params.addClassDescription(
      "Computes the \"area\" or dimension - 1 \"volume\" of a given "
      "boundary or boundaries in your mesh.");
  return params;
}

InterfaceAreaPostprocessor::InterfaceAreaPostprocessor(
    const InputParameters &parameters)
    : InterfaceIntegralPostprocessor(parameters) {}

void InterfaceAreaPostprocessor::threadJoin(const UserObject &y) {
  const InterfaceAreaPostprocessor &pps =
      static_cast<const InterfaceAreaPostprocessor &>(y);
  _integral_value += pps._integral_value;
}

Real InterfaceAreaPostprocessor::computeQpIntegral() { return 1.0; }
