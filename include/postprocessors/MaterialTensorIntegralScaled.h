//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MaterialTensorIntegral.h"

/**
 * Computes a volume integral of a material tensor component * scaled by
 * the value of another postprocessor.
 */
template <bool is_ad>
class MaterialTensorIntegralScaledTempl
    : public MaterialTensorIntegralTempl<is_ad> {
public:
  static InputParameters validParams();

  MaterialTensorIntegralScaledTempl(const InputParameters &parameters);
  virtual Real getValue() override;

protected:
  const PostprocessorValue &_scaling_factor_PP;
};

typedef MaterialTensorIntegralScaledTempl<false> MaterialTensorIntegralScaled;
typedef MaterialTensorIntegralScaledTempl<true> ADMaterialTensorIntegralScaled;
