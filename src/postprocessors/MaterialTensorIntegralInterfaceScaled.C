//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MaterialTensorIntegralInterfaceScaled.h"
#include "RankTwoScalarTools.h"

#include "metaphysicl/raw_type.h"

registerMooseObject("DeerApp", MaterialTensorIntegralInterfaceScaled);
registerMooseObject("DeerApp", ADMaterialTensorIntegralInterfaceScaled);

template <bool is_ad>
InputParameters
MaterialTensorIntegralInterfaceScaledTempl<is_ad>::validParams() {
  InputParameters params = InterfaceIntegralPostprocessor::validParams();
  params.addClassDescription("Computes an interface integral of "
                             "a component of a material tensor as specified by "
                             "the user-supplied indices.");
  params.addRequiredParam<MaterialPropertyName>(
      "rank_two_tensor", "The rank two material tensor name");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i", "index_i >= 0 & index_i <= 2",
      "The index i of ij for the tensor to output (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j", "index_j >= 0 & index_j <= 2",
      "The index j of ij for the tensor to output (0, 1, 2)");
  params.addParam<PostprocessorName>(
      "scaling_factor_PP",
      "a postprocessor used as scaling factor for the integral");
  return params;
}

template <bool is_ad>
MaterialTensorIntegralInterfaceScaledTempl<
    is_ad>::MaterialTensorIntegralInterfaceScaledTempl(const InputParameters
                                                           &parameters)
    : InterfaceIntegralPostprocessor(parameters),
      _tensor(
          getGenericMaterialProperty<RankTwoTensor, is_ad>("rank_two_tensor")),
      _i(getParam<unsigned int>("index_i")),
      _j(getParam<unsigned int>("index_j")),
      _scaling_factor_PP(
          isParamValid("scaling_factor_PP")
              ? &getPostprocessorValueByName(
                    getParam<PostprocessorName>("scaling_factor_PP"))
              : nullptr) {}

template <bool is_ad>
Real MaterialTensorIntegralInterfaceScaledTempl<is_ad>::computeQpIntegral() {
  return RankTwoScalarTools::component(MetaPhysicL::raw_value(_tensor[_qp]), _i,
                                       _j);
}

template <bool is_ad>
Real MaterialTensorIntegralInterfaceScaledTempl<is_ad>::getValue() {
  _integral_value = InterfaceIntegralPostprocessor::getValue();
  if (_scaling_factor_PP)
    _integral_value /= *_scaling_factor_PP;

  return _integral_value;
}

template class MaterialTensorIntegralInterfaceScaledTempl<false>;
template class MaterialTensorIntegralInterfaceScaledTempl<true>;
