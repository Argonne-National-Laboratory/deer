#include "TensorRateMaterial.h"

registerMooseObject("DeerApp", TensorRateMaterial);

InputParameters TensorRateMaterial::validParams() {
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>(
      "rank_two_tensor", "the tensorial material property");

  return params;
}

TensorRateMaterial::TensorRateMaterial(const InputParameters &parameters)
    : Material(parameters),
      _tensor(getMaterialPropertyByName<RankTwoTensor>(
          getParam<MaterialPropertyName>("rank_two_tensor"))),
      _tensor_old(getMaterialPropertyOld<RankTwoTensor>(
          getParam<MaterialPropertyName>("rank_two_tensor"))),
      _tensor_rate(declareProperty<RankTwoTensor>(
          getParam<MaterialPropertyName>("rank_two_tensor") + "_rate")) {}

void TensorRateMaterial::computeQpProperties() {
  _tensor_rate[_qp] = (_tensor[_qp] - _tensor_old[_qp]) / _dt;
}
