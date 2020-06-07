//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMInterfaceKernelPK1.h"

registerMooseObject("DeerApp", CZMInterfaceKernelPK1);

InputParameters CZMInterfaceKernelPK1::validParams() {
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredParam<unsigned int>(
      "component", "the component of the "
                   "displacement vector this kernel is working on:"
                   " component == 0, ==> X"
                   " component == 1, ==> Y"
                   " component == 2, ==> Z");
  params.set<bool>("_use_undisplaced_reference_points") = true;

  params.addRequiredCoupledVar("displacements",
                               "the string containing displacement variables");

  params.addClassDescription(
      "Interface kernel for use with cohesive zone models (CZMs) that "
      "compute traction as a function of the displacement jump");

  return params;
}

CZMInterfaceKernelPK1::CZMInterfaceKernelPK1(const InputParameters &parameters)
    : InterfaceKernel(parameters),
      _component(getParam<unsigned int>("component")),
      _ndisp(coupledComponents("displacements")), _disp_var(_ndisp),
      _disp_neighbor_var(_ndisp), _vars(_ndisp),
      _PK1traction(getMaterialPropertyByName<RealVectorValue>("PK1traction")),
      _dPK1traction_djumpglobal(
          getMaterialPropertyByName<RankTwoTensor>("dPK1traction_djumpglobal")),
      _dPK1traction_dF(
          getMaterialPropertyByName<RankThreeTensor>("dPK1traction_dF")) {
  if (getParam<bool>("use_displaced_mesh") == true)
    mooseError(
        "CZMInterfaceKernelPK1 cannot be used with use_displaced_mesh = true");

  for (unsigned int i = 0; i < _ndisp; ++i) {
    _disp_var[i] = coupled("displacements", i);
    _disp_neighbor_var[i] = coupled("displacements", i);
    _vars[i] = getVar("displacements", i);
  }
}

Real CZMInterfaceKernelPK1::computeQpResidual(Moose::DGResidualType type) {

  Real r = _PK1traction[_qp](_component);

  switch (type) {
  // [test_slave-test_master]*T where T represents the traction.
  case Moose::Element:
    r *= -_test[_i][_qp];
    break;
  case Moose::Neighbor:
    r *= _test_neighbor[_i][_qp];
    break;
  }
  return r;
}

Real CZMInterfaceKernelPK1::computeQpJacobian(Moose::DGJacobianType type) {
  // retrieve the diagonal Jacobian coefficient dependning on the displacement
  // component (_component) this kernel is working on
  Real jacsd = _dPK1traction_djumpglobal[_qp](_component, _component);
  Real jac = 0;
  switch (type) {
  case Moose::ElementElement: // Residual_sign -1  ddeltaU_ddisp sign -1;
    jac += _test[_i][_qp] * jacsd * _vars[_component]->phiFace()[_j][_qp];
    jac -= _test[_i][_qp] * JacLD(_component, /*neighbor=*/false);
    break;
  case Moose::ElementNeighbor: // Residual_sign -1  ddeltaU_ddisp sign 1;
    jac -=
        _test[_i][_qp] * jacsd * _vars[_component]->phiFaceNeighbor()[_j][_qp];
    jac -= _test[_i][_qp] * JacLD(_component, /*neighbor=*/true);
    break;
  case Moose::NeighborElement: // Residual_sign 1  ddeltaU_ddisp sign -1;
    jac -=
        _test_neighbor[_i][_qp] * jacsd * _vars[_component]->phiFace()[_j][_qp];
    jac += _test_neighbor[_i][_qp] * JacLD(_component, /*neighbor=*/false);
    break;
  case Moose::NeighborNeighbor: // Residual_sign 1  ddeltaU_ddisp sign 1;
    jac += _test_neighbor[_i][_qp] * jacsd *
           _vars[_component]->phiFaceNeighbor()[_j][_qp];
    jac += _test_neighbor[_i][_qp] * JacLD(_component, /*neighbor=*/true);
    break;
  }
  return jac;
}

Real CZMInterfaceKernelPK1::computeQpOffDiagJacobian(Moose::DGJacobianType type,
                                                     unsigned int jvar) {

  // find the displacement component associated to jvar
  unsigned int off_diag_component;
  for (off_diag_component = 0; off_diag_component < _ndisp;
       off_diag_component++)
    if (_disp_var[off_diag_component] == jvar)
      break;

  mooseAssert(off_diag_component < _ndisp,
              "CZMInterfaceKernelPK1::computeQpOffDiagJacobian wrong "
              "offdiagonal variable");

  Real jacsd = _dPK1traction_djumpglobal[_qp](_component, off_diag_component);
  Real jac = 0;

  switch (type) {
  case Moose::ElementElement: // Residual_sign -1  ddeltaU_ddisp sign -1;
    jac +=
        _test[_i][_qp] * jacsd * _vars[off_diag_component]->phiFace()[_j][_qp];
    jac -= _test[_i][_qp] * JacLD(off_diag_component, /*neighbor=*/false);
    break;
  case Moose::ElementNeighbor: // Residual_sign -1  ddeltaU_ddisp sign 1;
    jac -= _test[_i][_qp] * jacsd *
           _vars[off_diag_component]->phiFaceNeighbor()[_j][_qp];
    jac -= _test[_i][_qp] * JacLD(off_diag_component, /*neighbor=*/true);
    break;
  case Moose::NeighborElement: // Residual_sign 1  ddeltaU_ddisp sign -1;
    jac -= _test_neighbor[_i][_qp] * jacsd *
           _vars[off_diag_component]->phiFace()[_j][_qp];
    jac +=
        _test_neighbor[_i][_qp] * JacLD(off_diag_component, /*neighbor=*/false);
    break;
  case Moose::NeighborNeighbor: // Residual_sign 1  ddeltaU_ddisp sign 1;
    jac += _test_neighbor[_i][_qp] * jacsd *
           _vars[off_diag_component]->phiFaceNeighbor()[_j][_qp];
    jac +=
        _test_neighbor[_i][_qp] * JacLD(off_diag_component, /*neighbor=*/true);
    break;
  }
  return jac;
}

Real CZMInterfaceKernelPK1::JacLD(const unsigned int cc,
                                  const bool neighbor) const {
  Real jacld = 0;
  RealVectorValue phi;
  if (neighbor)
    phi = 0.5 * _vars[cc]->gradPhiFaceNeighbor()[_j][_qp];
  else
    phi = 0.5 * _vars[cc]->gradPhiFace()[_j][_qp];

  for (unsigned int j = 0; j < 3; j++)
    jacld += _dPK1traction_dF[_qp](_component, cc, j) * phi(j);
  return jacld;
}
