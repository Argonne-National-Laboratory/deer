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
#include "InterfaceIntegralPostprocessor.h"

/**
 * This postprocessor computes the cohesive zone area based on the provided
 * strain formulation
 */

class CZMAreaPostprocessor : public InterfaceIntegralPostprocessor
{
public:
  static InputParameters validParams();

  CZMAreaPostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  /// Base name of the material system
  const std::string _base_name;

  /// strain formulation
  enum class Strain
  {
    Small,
    Finite
  } _strain;

  /// the czm deformation gradient (if needed)
  const MaterialProperty<RankTwoTensor> * _F_czm;
};
