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

/**
 * MeshPropertiesMaterial uses the namespace EffectiveStressTools to compute
 * scalar values from Rank-2 tensors.
 */
class MeshPropertiesMaterial : public Material {
public:
  static InputParameters validParams();
  MeshPropertiesMaterial(const InputParameters &parameters);

protected:
  virtual void computeQpProperties() override;

  MaterialProperty<Real> &_h_min;
  MaterialProperty<Real> &_h_max;
  MaterialProperty<Real> &_jxw;
};
