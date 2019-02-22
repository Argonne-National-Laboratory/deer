#include "ComputeNEMLLargeStrain.h"

registerMooseObject("DeerApp", ComputeNEMLLargeStrain);

template <>
InputParameters
validParams<ComputeNEMLLargeStrain>()
{
  InputParameters params = validParams<ComputeNEMLStrainBase>();
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

ComputeNEMLLargeStrain::ComputeNEMLLargeStrain(const InputParameters &
                                               parameters)
  : ComputeNEMLStrainBase(parameters)
{

}

void
ComputeNEMLLargeStrain::computeQpStatefulProperties()
{
  RankTwoTensor F_np1 = RankTwoTensor(
        (*_grad_disp[0])[_qp],
        (*_grad_disp[1])[_qp],
        (*_grad_disp[2])[_qp]) + RankTwoTensor::Identity();

  RankTwoTensor F_n = RankTwoTensor(
        (*_grad_disp_old[0])[_qp],
        (*_grad_disp_old[1])[_qp],
        (*_grad_disp_old[2])[_qp]) + RankTwoTensor::Identity();

  RankTwoTensor Fi = F_np1.inverse();

  RankTwoTensor L = RankTwoTensor::Identity() - F_n * Fi;

  _strain_inc[_qp] = (L + L.transpose()) / 2.0;
  _mechanical_strain_inc[_qp] = (L + L.transpose()) / 2.0;
  _vorticity_inc[_qp] = (L - L.transpose()) / 2.0;
  
  RankFourTensor dL;
  dL.zero();
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        for (int m=0; m<3; m++) {
          for (int n=0; n<3; n++) {
            dL(i,j,m,n) += F_n(i,k) * Fi(k,m) * Fi(n,j);
          }
        }
      }
    }
  }

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        for (int l=0; l<3; l++) {
          _strain_grad[_qp](i,j,k,l) = 0.5 * (dL(i,j,k,l) + dL(j,i,k,l));
          _vorticity_grad[_qp](i,j,k,l) = 0.5 * (dL(i,j,k,l) - dL(j,i,k,l));          
        }
      }
    }
  }

  _ref_grad[_qp] = false;
}
