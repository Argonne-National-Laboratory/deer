#include "ComputeNEMLStrainBase.h"

template <>
InputParameters
validParams<ComputeNEMLStrainBase>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredCoupledVar("displacements",
                               "Displacement variables");
  params.addParam<bool>("large_kinematics", false, "Use large displacement kinematics.");
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

ComputeNEMLStrainBase::ComputeNEMLStrainBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _grad_disp(3),
    _strain_inc(declareProperty<RankTwoTensor>("strain_inc")),
    _mechanical_strain_inc(declareProperty<RankTwoTensor>("mechanical_strain_inc")),
    _vorticity_inc(declareProperty<RankTwoTensor>("vorticity_inc")),
    _def_grad(declareProperty<RankTwoTensor>("def_grad")),
    _def_grad_old(getMaterialPropertyOld<RankTwoTensor>("def_grad")),
    _df(declareProperty<RankTwoTensor>("df")),
    _ld(getParam<bool>("large_kinematics"))
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
  }

  // All others zero (so this will work naturally for plane strain problems)
  for (unsigned int i = _ndisp; i < 3; i++) {
    _disp[i] = &_zero;
    _grad_disp[i] = &_grad_zero;
  }
}

void
ComputeNEMLStrainBase::initQpStatefulProperties()
{
  _strain_inc[_qp].zero();
  _mechanical_strain_inc[_qp].zero();
  _vorticity_inc[_qp].zero();
  _def_grad[_qp] = RankTwoTensor::Identity();
  _df[_qp] = RankTwoTensor::Identity();
}

void ComputeNEMLStrainBase::computeProperties()
{
  if (_bnd) return;
  precalculate();
  DerivativeMaterialInterface<Material>::computeProperties();
}

void ComputeNEMLStrainBase::precalculate()
{

}

void ComputeNEMLStrainBase::computeQpProperties()
{
  _def_grad[_qp] = 
      (RankTwoTensor::Identity() + RankTwoTensor(
        (*_grad_disp[0])[_qp],
        (*_grad_disp[1])[_qp],
        (*_grad_disp[2])[_qp]));
}
