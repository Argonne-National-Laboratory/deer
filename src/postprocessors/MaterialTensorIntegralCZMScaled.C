//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MaterialTensorIntegralCZMScaled.h"
#include "RankTwoScalarTools.h"

#include "metaphysicl/raw_type.h"

registerMooseObject("DeerApp", MaterialTensorIntegralCZMScaled);
registerMooseObject("DeerApp", ADMaterialTensorIntegralCZMScaled);

template <bool is_ad>
InputParameters MaterialTensorIntegralCZMScaledTempl<is_ad>::validParams() {
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
      "A postprocessor used as scaling factor for the integral");
  params.addParam<bool>("normalize_integral_by_area", false,
                        "if true normalize the integral by the interface area. "
                        "This is done in addition to the scaling_factor_PP.");
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
  return params;
}

template <bool is_ad>
MaterialTensorIntegralCZMScaledTempl<is_ad>::
    MaterialTensorIntegralCZMScaledTempl(const InputParameters &parameters)
    : InterfaceIntegralPostprocessor(parameters),
      _tensor(
          getGenericMaterialProperty<RankTwoTensor, is_ad>("rank_two_tensor")),
      _i(getParam<unsigned int>("index_i")),
      _j(getParam<unsigned int>("index_j")),
      _scaling_factor_PP(
          isParamValid("scaling_factor_PP")
              ? &getPostprocessorValueByName(
                    getParam<PostprocessorName>("scaling_factor_PP"))
              : nullptr),
      _normalize_integral_by_area(getParam<bool>("normalize_integral_by_area")),
      _czm_area_mp(getMaterialPropertyByName<Real>("czm_area_mp")) {}

template <bool is_ad>
Real MaterialTensorIntegralCZMScaledTempl<is_ad>::computeQpIntegral() {
  return RankTwoScalarTools::component(MetaPhysicL::raw_value(_tensor[_qp]), _i,
                                       _j);
}

template <bool is_ad>
Real MaterialTensorIntegralCZMScaledTempl<is_ad>::computeIntegral() {
  Real sum = 0;
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    sum += _czm_area_mp[_qp] * _coord[_qp] * computeQpIntegral();
  return sum;
}

template <bool is_ad>
void MaterialTensorIntegralCZMScaledTempl<is_ad>::initialize() {
  InterfaceIntegralPostprocessor::initialize();
  _czm_area = 0;
}

template <bool is_ad>
void MaterialTensorIntegralCZMScaledTempl<is_ad>::execute() {
  InterfaceIntegralPostprocessor::execute();
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    _czm_area += _czm_area_mp[_qp];
}

template <bool is_ad>
void MaterialTensorIntegralCZMScaledTempl<is_ad>::finalize() {
  InterfaceIntegralPostprocessor::finalize();
  gatherSum(_czm_area);
  if (_scaling_factor_PP)
    _integral_value /= *_scaling_factor_PP;
  if (_normalize_integral_by_area)
    _integral_value /= _czm_area;
}

template class MaterialTensorIntegralCZMScaledTempl<false>;
template class MaterialTensorIntegralCZMScaledTempl<true>;
