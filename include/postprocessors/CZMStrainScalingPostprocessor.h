//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

class CZMStrainScalingPostprocessor;

template <> InputParameters validParams<CZMStrainScalingPostprocessor>();

/**
 * Computes teh ratio between the initial interface volume and initial interface
 * area.
 */
class CZMStrainScalingPostprocessor : public GeneralPostprocessor {
public:
  static InputParameters validParams();

  CZMStrainScalingPostprocessor(const InputParameters &parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual PostprocessorValue getValue() override;

protected:
  const PostprocessorValue &_a0;
  const PostprocessorValue &_v0;
};
