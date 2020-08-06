#include "TensorRateMaterial.h"

registerMooseObject("DeerApp", TensorRateMaterial);

InputParameters TensorRateMaterial::validParams() {
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>(
      "rank_two_tensor", "The tensor material property name");
  params.addParam<bool>(
      "convert_log_strain_to_eng_strain", false,
      "If true converts logarithmic strains to engineering strains");
  return params;
}

TensorRateMaterial::TensorRateMaterial(const InputParameters &parameters)
    : Material(parameters),
      _tensor(getMaterialPropertyByName<RankTwoTensor>(
          getParam<MaterialPropertyName>("rank_two_tensor"))),
      _tensor_old(getMaterialPropertyOld<RankTwoTensor>(
          getParam<MaterialPropertyName>("rank_two_tensor"))),
      _tensor_rate(declareProperty<RankTwoTensor>(
          getParam<MaterialPropertyName>("rank_two_tensor") + "_rate")),
      _log(getParam<bool>("convert_log_strain_to_eng_strain")) {}

void TensorRateMaterial::computeQpProperties() {
  RankTwoTensor t = _tensor[_qp];
  RankTwoTensor t_old = _tensor_old[_qp];
  if (_log)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        t(i, j) = std::exp(t(i, j)) - 1;
        t_old(i, j) = std::exp(t_old(i, j)) - 1;
      }
  _tensor_rate[_qp] = (t - t_old) / _dt;
}
