//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMStrainScalingPostprocessor.h"

registerMooseObject("DeerApp", CZMStrainScalingPostprocessor);

defineLegacyParams(CZMStrainScalingPostprocessor);

InputParameters CZMStrainScalingPostprocessor::validParams() {
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<PostprocessorName>("A0", "Initial Area");
  params.addRequiredParam<PostprocessorName>("V0", "Volume");
  params.addClassDescription("Computes the sum of two postprocessors");
  // params.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;
  return params;
}

CZMStrainScalingPostprocessor::CZMStrainScalingPostprocessor(
    const InputParameters &parameters)
    : GeneralPostprocessor(parameters), _a0(getPostprocessorValue("A0")),
      _v0(getPostprocessorValue("V0")) {}

void CZMStrainScalingPostprocessor::initialize() {}

void CZMStrainScalingPostprocessor::execute() {}

PostprocessorValue CZMStrainScalingPostprocessor::getValue() {
  return _v0 / _a0;
}
