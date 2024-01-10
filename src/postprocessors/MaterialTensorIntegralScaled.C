//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MaterialTensorIntegralScaled.h"

registerMooseObject("DeerApp", MaterialTensorIntegralScaled);
registerMooseObject("DeerApp", ADMaterialTensorIntegralScaled);

template <bool is_ad>
InputParameters
MaterialTensorIntegralScaledTempl<is_ad>::validParams()
{
  InputParameters params = MaterialTensorIntegralTempl<is_ad>::validParams();
  params.addRequiredParam<PostprocessorName>(
      "scaling_factor_PP", "A postprocessor used as scaling factor for the integral");
  params.addClassDescription("Computes a volume integral of a material tensor component scaled "
                             "by the value of another postprocessor.");
  return params;
}

template <bool is_ad>
MaterialTensorIntegralScaledTempl<is_ad>::MaterialTensorIntegralScaledTempl(
    const InputParameters & parameters)
  : MaterialTensorIntegralTempl<is_ad>(parameters),
    _scaling_factor_PP(this->getPostprocessorValueByName(
        this->template getParam<PostprocessorName>("scaling_factor_PP")))
{
}

template <bool is_ad>
PostprocessorValue
MaterialTensorIntegralScaledTempl<is_ad>::getValue()
{
  this->_integral_value = MaterialTensorIntegralTempl<is_ad>::getValue();
  this->_integral_value /= _scaling_factor_PP;
  return this->_integral_value;
}

template class MaterialTensorIntegralScaledTempl<false>;
template class MaterialTensorIntegralScaledTempl<true>;
