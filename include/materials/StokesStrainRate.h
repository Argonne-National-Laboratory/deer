//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

/**
 * Calculate the simple small strain rate
 */
class StokesStrainRate : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  StokesStrainRate(const InputParameters & parameters);

protected:
  /// @brief  Calculate the strain
  void computeQpProperties() override;

  /// The strain rate
  ADMaterialProperty<RankTwoTensor> & _strain_rate;

  /// The coupled velocity gradient
  const ADVectorVariableGradient & _grad_vel;
};