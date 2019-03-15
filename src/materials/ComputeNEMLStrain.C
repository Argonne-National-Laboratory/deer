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

  RankTwoTensor F_n = _def_grad_inv_old[_qp].inverse();
  
  RankTwoTensor L = RankTwoTensor::Identity() - F_n * _def_grad_inv[_qp];

  _strain_inc[_qp] = (L + L.transpose()) / 2.0;
  _mechanical_strain_inc[_qp] = (L + L.transpose()) / 2.0;
  _vorticity_inc[_qp] = (L - L.transpose()) / 2.0;
}
