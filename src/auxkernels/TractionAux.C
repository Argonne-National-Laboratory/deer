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

MooseEnum traction_component() {
  return MooseEnum("normal shear1 shear2 shear_norm");
}

template <> InputParameters validParams<TractionAux>() {
  InputParameters params = validParams<MaterialAuxBase<>>();
  params.addRequiredParam<MooseEnum>("scalar_type", traction_component(),
                                     "Type of scalar output");

  params.addClassDescription("Compute traction components on a boundary and "
                             "save it to an aux variable");

  return params;
}

TractionAux::TractionAux(const InputParameters &parameters)
    : MaterialAuxBase<RankTwoTensor>(parameters), _normals(_assembly.normals()),
      _scalar_type(getParam<MooseEnum>("scalar_type")) {
  if (!this->boundaryRestricted())
    mooseError("TractionAux Kernle must be boundary restricted ");
}

Real TractionAux::getRealValue() {

  RealVectorValue local_normal(1, 0, 0);

  // compute rotation matrix required to aligne the normal with the 1,0,0 vector
  RealTensorValue RotationGlobal2Local =
      RotationMatrix::rotVec1ToVec2(_normals[_qp], local_normal);

  // rotate the stress tensor
  RankTwoTensor stress_R = _prop[_qp];
  stress_R.rotate(RotationGlobal2Local);

  // compute the traction
  RealVectorValue traction = stress_R * local_normal;

  // extract components
  Real val;
  switch (_scalar_type) {
  case 0:
    val = traction(0);
    break;
  case 1:
    val = traction(1);
    break;
  case 2:
    val = traction(2);
    break;
  case 3:
    val = std::sqrt(traction(1) * traction(1) + traction(2) * traction(2));
    break;
  }

  return val;
}
