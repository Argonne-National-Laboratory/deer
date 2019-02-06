#include "ComputeNEMLSmallStrain.h"

registerMooseObject("DeerApp", ComputeNEMLSmallStrain);

template <>
InputParameters
validParams<ComputeNEMLSmallStrain>()
{
  InputParameters params = validParams<ComputeNEMLStrainBase>();
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

ComputeNEMLSmallStrain::ComputeNEMLSmallStrain(const InputParameters &
                                               parameters)
  : ComputeNEMLStrainBase(parameters)
{

}

void
ComputeNEMLSmallStrain::computeProperties()
{
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    RankTwoTensor grad_inc(
        (*_grad_disp[0])[_qp] - (*_grad_disp_old[0])[_qp],
        (*_grad_disp[1])[_qp] - (*_grad_disp_old[1])[_qp],
        (*_grad_disp[2])[_qp] - (*_grad_disp_old[2])[_qp]);
    _strain_inc[_qp] = (grad_inc + grad_inc.transpose()) / 2.0;
    _mechanical_strain_inc[_qp] = (grad_inc + grad_inc.transpose()) / 2.0;

    _vorticity_inc[_qp].zero();

    _strain_grad[_qp] =
        RankFourTensor(RankFourTensor::initIdentitySymmetricFour);
    _vorticity_grad[_qp].zero();
  }
}
