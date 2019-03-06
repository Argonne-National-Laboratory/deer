#include "ComputeNEMLLargeStrain.h"

registerMooseObject("DeerApp", ComputeNEMLLargeStrain);

template <>
InputParameters
validParams<ComputeNEMLLargeStrain>()
{
  InputParameters params = validParams<ComputeNEMLStrainBase>();
  return params;
}

ComputeNEMLLargeStrain::ComputeNEMLLargeStrain(const InputParameters &
                                               parameters)
  : ComputeNEMLStrainBase(parameters)
{

}

void
ComputeNEMLLargeStrain::computeQpProperties()
{
  ComputeNEMLStrainBase::computeQpProperties();

  RankTwoTensor F_n = _def_grad_inv_old[_qp].inverse();
  
  RankTwoTensor L = RankTwoTensor::Identity() - F_n * _def_grad_inv[_qp];

  _strain_inc[_qp] = (L + L.transpose()) / 2.0;
  _mechanical_strain_inc[_qp] = (L + L.transpose()) / 2.0;
  _vorticity_inc[_qp] = (L - L.transpose()) / 2.0;

  RankFourTensor dL = F_n.mixedProductIkJl(RankTwoTensor::Identity());

  // MOOSE doesn't have the skew-symmetric identity baked in
  RankFourTensor dD;
  RankFourTensor dW;
  dD.zero();
  dW.zero();
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        for (int l=0; l<3; l++) {
          dD(i,j,k,l) = 0.5*((i==l)*(j==k) + (i==k)*(j==l)); // this is just initIdentitySymmetricFour
          dW(i,j,k,l) = 0.5*((i==l)*(j==k) - (i==k)*(j==l));
        }
      }
    }
  }

  _strain_grad[_qp] = dD * dL;
  _vorticity_grad[_qp] = dW * dL;
}
