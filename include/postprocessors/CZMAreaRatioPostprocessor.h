//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "CZMAreaPostprocessor.h"

/**
 * This postprocessor computes the cohesive zon area ratio, i.e. A/A0
 */

class CZMAreaRatioPostprocessor : public CZMAreaPostprocessor {
public:
  static InputParameters validParams();

  CZMAreaRatioPostprocessor(const InputParameters &parameters);

protected:
  virtual Real getValue() override;
};
