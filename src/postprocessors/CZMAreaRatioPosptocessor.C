//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMAreaRatioPostprocessor.h"

registerMooseObject("DeerApp", CZMAreaRatioPostprocessor);

InputParameters
CZMAreaRatioPostprocessor::validParams()
{
  InputParameters params = CZMAreaPostprocessor::validParams();
  params.addClassDescription("Computes the area ratio, i.e. A/A0, of a cohesive zone.");
  return params;
}

CZMAreaRatioPostprocessor::CZMAreaRatioPostprocessor(const InputParameters & parameters)
  : CZMAreaPostprocessor(parameters)
{
}

PostprocessorValue
CZMAreaRatioPostprocessor::getValue() const
{
  return CZMAreaPostprocessor::getValue() / _interface_primary_area;
}
