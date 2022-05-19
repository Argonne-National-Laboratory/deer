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
class MaterialTensorIntegralInterfaceScaledTempl : public InterfaceIntegralPostprocessor
{
public:
  static InputParameters validParams();

  MaterialTensorIntegralInterfaceScaledTempl(const InputParameters & parameters);
  virtual Real getValue() override;

protected:
  virtual Real computeQpIntegral();

private:
  const GenericMaterialProperty<RankTwoTensor, is_ad> & _tensor;
  const unsigned int _i;
  const unsigned int _j;
  const PostprocessorValue * _scaling_factor_PP;
};

typedef MaterialTensorIntegralInterfaceScaledTempl<false> MaterialTensorIntegralInterfaceScaled;
typedef MaterialTensorIntegralInterfaceScaledTempl<true> ADMaterialTensorIntegralInterfaceScaled;
