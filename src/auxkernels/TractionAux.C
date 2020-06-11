//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Assembly.h"
#include "MooseMesh.h"
#include "RotationMatrix.h"
#include "TractionAux.h"

#include "libmesh/dense_matrix.h"

registerMooseObject("DeerApp", TractionAux);

MooseEnum scalar_type("normal shear1 shear2 shear_norm X Y Z");

InputParameters TractionAux::validParams() {
  InputParameters params = MaterialAuxBase<RankTwoTensor>::validParams();
  params.addRequiredParam<MooseEnum>(
      "scalar_type", scalar_type,
      "Type of scalar output, possible options are: normal, shear1, shear2, "
      "shear_norm, X, Y, Z");
  params.addParam<bool>("PK1", false, "traction value on the undeformed area.");
  params.addParam<bool>("large_kinematics", true,
                        "use the updated geometry to compute the traction");
  params.addClassDescription("Compute a traction component on a boundary and "
                             "save it to an aux variable eihter in global or "
                             "interface natural cooridantes");
  params.addRequiredCoupledVar(
      "displacements",
      "The string of displacements suitable for the problem statement");
  return params;
}

TractionAux::TractionAux(const InputParameters &parameters)
    : MaterialAuxBase<RankTwoTensor>(parameters), _normals(_assembly.normals()),
      _PK1(getParam<bool>("PK1")),
      _large_kinematics(getParam<bool>("large_kinematics")),
      _ndisp(coupledComponents("displacements")),
      _scalar_type(getParam<MooseEnum>("scalar_type").getEnum<scalarType>()) {
  if (!this->boundaryRestricted())
    mooseError("TractionAux Kernel must be boundary restricted!");
}

void TractionAux::initialSetup() {
  for (unsigned int i = 0; i < _ndisp; ++i)
    _grad_disp.push_back(&coupledGradient("displacements", i));

  // All others zero (so this will work naturally for 2D and 1D problems)
  for (unsigned int i = _ndisp; i < 3; i++)
    _grad_disp.push_back(&_grad_zero);
}

Real TractionAux::getRealValue() {

  RealVectorValue local_normal(1, 0, 0);
  RealVectorValue normal = _normals[_qp];
  Real dadA = 1.;
  if (_large_kinematics) {
    RankTwoTensor F =
        (RankTwoTensor::Identity() + RankTwoTensor((*_grad_disp[0])[_qp],
                                                   (*_grad_disp[1])[_qp],
                                                   (*_grad_disp[2])[_qp]));
    RankTwoTensor R;
    F.getRUDecompositionRotation(R);
    normal = R * _normals[_qp];
    if (_PK1)
      dadA = F.det() * (F.inverse().transpose() * _normals[_qp]).norm();
  }

  // compute rotation matrix required to align the normal with the {1,0,0}
  // vector
  RealTensorValue RotationGlobal2Local =
      RotationMatrix::rotVec1ToVec2(normal, local_normal);

  // rotate the stress tensor
  RankTwoTensor stress_R = _prop[_qp];
  stress_R.rotate(RotationGlobal2Local);

  // compute the traction
  RealVectorValue traction = stress_R * local_normal;
  // scale traction on original area
  traction *= dadA;

  // extract components
  Real val = 0;
  switch (_scalar_type) {
  case scalarType::normal:
    val = traction(0);
    break;
  case scalarType::shear1:
    val = traction(1);
    break;
  case scalarType::shear2:
    val = traction(2);
    break;
  case scalarType::shear_norm:
    val = std::sqrt(traction(1) * traction(1) + traction(2) * traction(2));
    break;
  case scalarType::X:
    traction = RotationGlobal2Local.transpose() * traction;
    val = traction(0);
    break;
  case scalarType::Y:
    traction = RotationGlobal2Local.transpose() * traction;
    val = traction(1);
    break;
  case scalarType::Z:
    traction = RotationGlobal2Local.transpose() * traction;
    val = traction(2);
    break;
  }

  return val;
}
