#include "ComputeNEMLStrainBase.h"

template <>
InputParameters
validParams<ComputeNEMLStrainBase>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredCoupledVar("displacements",
                               "Displacement variables");

  return params;
}

ComputeNEMLStrainBase::ComputeNEMLStrainBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _grad_disp(3),
    _disp_old(3),
    _grad_disp_old(3),
    _strain_inc(declareProperty<RankTwoTensor>("strain_inc")),
    _mechanical_strain_inc(declareProperty<RankTwoTensor>("mechanical_strain_inc")),
    _vorticity_inc(declareProperty<RankTwoTensor>("vorticity_inc")),
    _strain_grad(declareProperty<RankFourTensor>("strain_grad")),
    _vorticity_grad(declareProperty<RankFourTensor>("vorticity_grad")),
    _ref_grad(declareProperty<bool>("strain_reference_grad"))
{

}

void
ComputeNEMLStrainBase::initialSetup()
{
  // Enforce consistency
  if (_ndisp != _mesh.dimension()) {
    paramError("displacements", 
               "Number of displacements must match problem dimension.");
  }

  for (unsigned int i = 0; i < _ndisp; i++) {
    _disp[i] = &coupledValue("displacements", i);
    _grad_disp[i] = &coupledGradient("displacements", i);
    _disp_old[i] = &coupledValueOld("displacements", i);
    _grad_disp_old[i] = &coupledGradientOld("displacements", i);
  }

  // All others zero (so this will work naturally for plane strain problems)
  for (unsigned int i = _ndisp; i < 3; i++) {
    _disp[i] = &_zero;
    _grad_disp[i] = &_grad_zero;
    _disp_old[i] = &_zero;
    _grad_disp_old[i] = &_grad_zero;
  }
}

void
ComputeNEMLStrainBase::initQpStatefulProperties()
{
  _strain_inc[_qp].zero();
  _mechanical_strain_inc[_qp].zero();
  _vorticity_inc[_qp].zero();
  _strain_grad[_qp].zero();
  _vorticity_grad[_qp].zero();
}

void ComputeNEMLStrainBase::computeProperties()
{
  precalculate();
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    computeQpStatefulProperties();
  }
}

void ComputeNEMLStrainBase::precalculate()
{

}
