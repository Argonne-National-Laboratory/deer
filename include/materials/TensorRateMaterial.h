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
 * Material class calcualting a tensor rate component by componet
 */
class TensorRateMaterial : public Material {
public:
  static InputParameters validParams();
  TensorRateMaterial(const InputParameters &parameters);

protected:
  virtual void computeQpProperties() override;
  const MaterialProperty<RankTwoTensor> &_tensor;
  const MaterialProperty<RankTwoTensor> &_tensor_old;
  MaterialProperty<RankTwoTensor> &_tensor_rate;
};
