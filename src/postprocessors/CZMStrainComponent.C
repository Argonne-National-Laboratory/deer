//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMStrainComponent.h"

registerMooseObject("DeerApp", CZMStrainComponent);

InputParameters
CZMStrainComponent::validParams()
{
  InputParameters params = CZMAreaRatioPostprocessor::validParams();
  params.addClassDescription(
      "Postprocessor computing the cohesive zon strain contrinution to an RVE");
  params.addRequiredParam<MaterialPropertyName>("rank_two_tensor",
                                                "The rank two material tensor name");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i",
      "index_i >= 0 & index_i <= 2",
      "The index i of ij for the tensor to output (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j",
      "index_j >= 0 & index_j <= 2",
      "The index j of ij for the tensor to output (0, 1, 2)");
  params.addRequiredParam<PostprocessorName>(
      "initial_bulk_volume_pp", "A postprocessor used as scaling factor for the integral");
  return params;
}

CZMStrainComponent::CZMStrainComponent(const InputParameters & parameters)
  : CZMAreaRatioPostprocessor(parameters),
    _tensor(getMaterialPropertyByName<RankTwoTensor>(
        _base_name + getParam<MaterialPropertyName>("rank_two_tensor"))),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j")),
    _initial_bulk_volume_pp(
        getPostprocessorValueByName(getParam<PostprocessorName>("initial_bulk_volume_pp")))
{
}

void
CZMStrainComponent::initialize()
{
  CZMAreaRatioPostprocessor::initialize();
  _normalized_strain_component = 0;
}

PostprocessorValue
CZMStrainComponent::getValue() const
{
  return _normalized_strain_component / (_initial_bulk_volume_pp * CZMAreaRatioPostprocessor::getValue());
}

void
CZMStrainComponent::finalize()
{
  gatherSum(_normalized_strain_component);
}

void
CZMStrainComponent::execute()
{
  CZMAreaRatioPostprocessor::execute();
  _normalized_strain_component += computeStrainIntegral();
}

Real
CZMStrainComponent::computeStrainIntegral()
{
  Real sum = 0;
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    sum += _JxW[_qp] * _coord[_qp] * _tensor[_qp](_i, _j);
  return sum;
}

void
CZMStrainComponent::threadJoin(const UserObject & y)
{
  InterfaceIntegralPostprocessor::threadJoin(y);
  const CZMStrainComponent & pps = static_cast<const CZMStrainComponent &>(y);
  _normalized_strain_component += pps._normalized_strain_component;
}
