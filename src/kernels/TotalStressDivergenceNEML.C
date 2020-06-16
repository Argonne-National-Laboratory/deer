#include "TotalStressDivergenceNEML.h"

registerMooseObject("DeerApp", TotalStressDivergenceNEML);

InputParameters TotalStressDivergenceNEML::validParams() {
  InputParameters params = Kernel::validParams();

  params.addRequiredParam<unsigned int>("component",
                                        "Which direction this kernel acts in");
  params.addRequiredCoupledVar("displacements", "The displacement components");
  params.addParam<bool>("large_kinematics", false,
                        "Use large displacement kinematics.");
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

TotalStressDivergenceNEML::TotalStressDivergenceNEML(const InputParameters &parameters)
    : DerivativeMaterialInterface<Kernel>(parameters),
      _ld(getParam<bool>("large_kinematics")),
      _component(getParam<unsigned int>("component")),
      _ndisp(coupledComponents("displacements")), _disp_nums(_ndisp),
      _disp_vars(_ndisp), _grad_disp(_ndisp),
      _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
      _material_jacobian(
          getMaterialPropertyByName<RankFourTensor>("material_jacobian")),
      _inv_def_grad(getMaterialPropertyByName<RankTwoTensor>("inv_def_grad")),
      _detJ(getMaterialPropertyByName<Real>("detJ")),
      _df(getMaterialPropertyByName<RankTwoTensor>("df"))
{}

void 
TotalStressDivergenceNEML::initialSetup() {
  for (unsigned int i = 0; i < _ndisp; i++) {
    _disp_nums[i] = coupled("displacements", i);
    _disp_vars[i] = getVar("displacements", i);
    _grad_disp[i] = &coupledGradient("displacements", i);
  }
}

Real TotalStressDivergenceNEML::computeQpResidual() {
  if (_ld) {
    return largeDeformationResidual(_component, _grad_test[_i][_qp]);
  }
  else {
    return smallDeformationResidual(_component, _grad_test[_i][_qp]);
  }
}

Real
TotalStressDivergenceNEML::largeDeformationResidual(unsigned int i,
                                                    const RealGradient & grad_phi)
{
  return _detJ[_qp] * _stress[_qp].row(i) * 
      (_inv_def_grad[_qp].transpose() * grad_phi);
}

Real
TotalStressDivergenceNEML::smallDeformationResidual(unsigned int i,
                                                    const RealGradient & grad_phi)
{
  return _stress[_qp].row(i) * grad_phi;
}

Real TotalStressDivergenceNEML::computeQpJacobian() {
  Real value = 0.0;
  
  if (_ld) {
    value += largeDeformationMatJac(_component, _component,
                                    _grad_test[_i][_qp],
                                    _grad_phi[_j][_qp]);
    value += largeDeformationGeoJac(_component, _component,
                                    _grad_test[_i][_qp],
                                    _grad_phi[_j][_qp]);
  }
  else {
    value += smallDeformationMatJac(_component, _component,
                                    _grad_test[_i][_qp],
                                    _grad_phi[_j][_qp]);
  }

  return value;
}

Real TotalStressDivergenceNEML::computeQpOffDiagJacobian(unsigned int jvar) {
  Real value = 0.0;

  for (unsigned int cc = 0; cc < _ndisp; cc++) {
    if (jvar == _disp_nums[cc]) {
      if (_ld) {
        value += largeDeformationMatJac(_component, cc,
                                        _grad_test[_i][_qp],
                                        _disp_vars[cc]->gradPhi()[_j][_qp]);
        value += largeDeformationGeoJac(_component, _component,
                                        _grad_test[_i][_qp],
                                        _disp_vars[cc]->gradPhi()[_j][_qp]);
      }
      else {
        value += smallDeformationMatJac(_component, cc,
                                        _grad_test[_i][_qp],
                                        _disp_vars[cc]->gradPhi()[_j][_qp]);
      }
      break;
    }
  }

  return value;
}

Real 
TotalStressDivergenceNEML::smallDeformationMatJac(
    unsigned int i, unsigned int k, const RealGradient & grad_phi,
    const RealGradient & grad_psi)
{
  Real value = 0.0;
  for (unsigned int j = 0; j < _ndisp; j++) {
    for (unsigned int l = 0; l < _ndisp; l++) {
      value += _material_jacobian[_qp](i,j,k,l) * grad_phi(j) * grad_psi(l);
    }
  }
  return value;
}

Real
TotalStressDivergenceNEML::largeDeformationMatJac(
    unsigned int i, unsigned int k, const RealGradient & grad_phi,
    const RealGradient & grad_psi)
{
  auto GPhi = _inv_def_grad[_qp].transpose() * grad_phi;
  auto GPsi = _inv_def_grad[_qp].transpose() * grad_psi;

  Real value = 0.0;

  for (unsigned int j = 0; j < _ndisp; j++) {
    for (unsigned int m = 0; m < _ndisp; m++) {
      for (unsigned int n = 0; n < _ndisp; n++) {
        value += _detJ[_qp] * _material_jacobian[_qp](i,j,m,n) * 
            _df[_qp](m,k) * GPsi(n) * GPhi(j);
      }
    }
  }

  return value;
}

Real
TotalStressDivergenceNEML::largeDeformationGeoJac(
    unsigned int i, unsigned int k, const RealGradient & grad_phi,
    const RealGradient & grad_psi)
{
  auto GPhi = _inv_def_grad[_qp].transpose() * grad_phi;
  auto GPsi = _inv_def_grad[_qp].transpose() * grad_psi;
  
  Real value = 0.0;

  for (unsigned int j = 0; j < _ndisp; j++) {
    value += _detJ[_qp] * _stress[_qp](i,j) * 
        (GPsi(j) * GPhi(k) - GPsi(k) * GPhi(j));
  }

  return value;
}
