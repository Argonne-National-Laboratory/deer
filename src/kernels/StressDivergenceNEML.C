#include "StressDivergenceNEML.h"

registerMooseObject("DeerApp", StressDivergenceNEML);

template <>
InputParameters
validParams<StressDivergenceNEML>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredParam<unsigned int>("component", 
                                        "Which direction this kernel acts in");
  params.addRequiredCoupledVar("displacements",
                               "The displacement components");
  return params;
}

StressDivergenceNEML::StressDivergenceNEML(const InputParameters & parameters) 
  : DerivativeMaterialInterface<Kernel>(parameters),
    _ld(getParam<bool>("use_displaced_mesh")),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_nums(_ndisp),
    _disp_vars(_ndisp),
    _grad_disp(_ndisp),
    _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
    _material_strain_jacobian(getMaterialPropertyByName<RankFourTensor>("material_strain_jacobian")),
    _material_vorticity_jacobian(getMaterialPropertyByName<RankFourTensor>("material_vorticity_jacobian")),
    _strain_grad(getMaterialPropertyByName<RankFourTensor>("strain_grad")),
    _vorticity_grad(getMaterialPropertyByName<RankFourTensor>("vorticity_grad")),
    _inv_def_grad(getMaterialPropertyByName<RankTwoTensor>("def_grad_inv")),
    _inv_def_grad_old(getMaterialPropertyOld<RankTwoTensor>("def_grad_inv"))
{

}

void StressDivergenceNEML::initialSetup()
{
  for (unsigned int i = 0; i < _ndisp; i++) {
    _disp_nums[i] = coupled("displacements", i);
    _disp_vars[i] = getVar("displacements", i);
    _grad_disp[i] = &coupledGradient("displacements", i);
  }
}

void StressDivergenceNEML::precalculateResidual()
{

}

void StressDivergenceNEML::precalculateJacobian()
{

}

void StressDivergenceNEML::precalculateOffDiagJacobian(unsigned int jvar)
{

}

Real StressDivergenceNEML::computeQpResidual()
{
  return _stress[_qp].row(_component) * _grad_test[_i][_qp];
}

Real StressDivergenceNEML::computeQpJacobian()
{
  Real value = 0.0;
  
  if (not _ld) {
    value += matJacobianComponent(_material_strain_jacobian[_qp],
                                  _material_vorticity_jacobian[_qp],
                                  _strain_grad[_qp],
                                  _vorticity_grad[_qp],
                                  _component, _component,
                                  _grad_test[_i][_qp], _grad_phi[_j][_qp]);
  }
  else {
    value += geomJacobianComponent(_component, _component,
                                   _grad_test[_i][_qp], _grad_phi[_j][_qp],
                                   _stress[_qp]);
  }
  
  return value;
}

Real StressDivergenceNEML::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real value = 0.0;

  for (unsigned int cc = 0; cc < _ndisp; cc++) {
    if (jvar == _disp_nums[cc]) {
      if (not _ld) {
        value += matJacobianComponent(_material_strain_jacobian[_qp],
                                      _material_vorticity_jacobian[_qp],
                                      _strain_grad[_qp],
                                      _vorticity_grad[_qp],
                                      _component, cc,
                                      _grad_test[_i][_qp], 
                                      _disp_vars[cc]->gradPhi()[_j][_qp]);
      }
      else {
        value += geomJacobianComponent(_component, cc,
                                       _grad_test[_i][_qp], 
                                       _disp_vars[cc]->gradPhi()[_j][_qp],
                                       _stress[_qp]);
      }
    }
  }

  return value;
}

Real StressDivergenceNEML::matJacobianComponent(
    const RankFourTensor & A, const RankFourTensor & B, 
    const RankFourTensor & E, const RankFourTensor & W,
    unsigned int i, unsigned int m,
    const RealGradient & grad_psi, const RealGradient & grad_phi)
{
  // In index notation this is: 
  // grad_psi(i,j) * (A(i,j,k,l) * E(k,l,m,n) + B(i,j,k,l) * W(k,l,m,n)) *
  // grad_phi(m,n)
  
  // For now we will just do it that way, need to work out a better, more 
  // vectorized method to compute it

  Real value = 0.0;

  for (int j = 0; j < _mesh.dimension(); j++) {
    for (int k = 0; k < _mesh.dimension(); k++) {
      for (int l = 0; l < _mesh.dimension(); l++) {
        for (int n = 0; n < _mesh.dimension(); n++) {
          value += grad_psi(j) * (A(i,j,k,l) * E(k,l,m,n) 
                               + B(i,j,k,l) * W(k,l,m,n)) * grad_phi(n);
        }
      }
    }
  }

  return value;
}

Real StressDivergenceNEML::geomJacobianComponent(
    unsigned int i, unsigned int m,
    const RealGradient & grad_psi, const RealGradient & grad_phi,
    const RankTwoTensor & stress)
{
  RankTwoTensor F_n = _inv_def_grad_old[_qp].inverse();
  RankTwoTensor prod = F_n * _inv_def_grad[_qp];

  Real value = 0.0;

  for (int j=0; j<3; j++) {
    value += grad_phi(m) * stress(i,j) * grad_psi(j);
    value -= grad_phi(j) * stress(i,j) * grad_psi(m);
    for (int k=0; k<3; k++) {
      for (int n=0; n<3; n++) {
        value += (_material_strain_jacobian[_qp](i,j,k,n) 
                  + _material_vorticity_jacobian[_qp](i,j,k,n)) * prod(k,m) * 
            grad_phi(n) * grad_psi(j);
      }
    }
  }

  return value;
}
