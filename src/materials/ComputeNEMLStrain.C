#include "ComputeNEMLStrain.h"

registerMooseObject("DeerApp", ComputeNEMLStrain);

template <>
InputParameters
validParams<ComputeNEMLStrain>()
{
  InputParameters params = validParams<ComputeNEMLStrainBase>();
  return params;
}

ComputeNEMLStrain::ComputeNEMLStrain(const InputParameters &
                                               parameters)
  : ComputeNEMLStrainBase(parameters)
{

}

void
ComputeNEMLStrain::computeQpProperties()
{
  ComputeNEMLStrainBase::computeQpProperties();

  RankTwoTensor L;
  
  if (_ld) {
    L = RankTwoTensor::Identity() - _def_grad_old[_qp] * _def_grad[_qp].inverse();
    _df[_qp] = -L + RankTwoTensor::Identity();
  }
  else {
    L = _def_grad[_qp] - _def_grad_old[_qp];
    _df[_qp] = RankTwoTensor::Identity();
  }

  _strain_inc[_qp] = (L + L.transpose()) / 2.0;
  _mechanical_strain_inc[_qp] = (L + L.transpose()) / 2.0 - eigenstrainIncrement();
  _vorticity_inc[_qp] = (L - L.transpose()) / 2.0;
}
