//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceIntegralPostprocessor.h"
#include "RankTwoTensor.h"

/**
 *  This postprocessor computes an interface integral of a component of a
 * material tensor as specified by the user-supplied indices. Optionally the
 * resulting integral can be scaled by the value of another postprocessor.
 */
template <bool is_ad>
class MaterialTensorIntegralCZMScaledTempl
    : public InterfaceIntegralPostprocessor {
public:
  static InputParameters validParams();

  MaterialTensorIntegralCZMScaledTempl(const InputParameters &parameters);
  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;

protected:
  Real computeIntegral() override;
  Real computeQpIntegral() override;

private:
  const GenericMaterialProperty<RankTwoTensor, is_ad> &_tensor;
  const unsigned int _i;
  const unsigned int _j;
  const PostprocessorValue *_scaling_factor_PP;
  const bool _normalize_integral_by_area;
  Real _czm_area = 0;
  const MaterialProperty<Real> &_czm_area_mp;
};

typedef MaterialTensorIntegralCZMScaledTempl<false>
    MaterialTensorIntegralCZMScaled;
typedef MaterialTensorIntegralCZMScaledTempl<true>
    ADMaterialTensorIntegralCZMScaled;
