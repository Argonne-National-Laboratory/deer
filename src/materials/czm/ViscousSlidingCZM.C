//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceValueTools.h"
#include "RotationMatrix.h"
#include "ViscousSlidingCZM.h"

registerMooseObject("DeerApp", ViscousSlidingCZM);

template <> InputParameters validParams<ViscousSlidingCZM>() {
  InputParameters params = validParams<PureElasticCZM>();
  params.addClassDescription(
      "time dependent interface damage model with linear interpolation");
  params.addRequiredParam<Real>(
      "shear_viscosity", "interface shear viscosity [Pressure time / length]");
  params.addRequiredParam<Real>("interface_thickness",
                                "the interface thickness [length]");
  return params;
}

ViscousSlidingCZM::ViscousSlidingCZM(const InputParameters &parameters)
    : PureElasticCZM(parameters),
      _shear_viscosity(getParam<Real>("shear_viscosity")),
      _traction_old(getMaterialPropertyOld<RealVectorValue>("traction")),
      _displacement_jump_dot(
          declareProperty<RealVectorValue>("displacement_jump_dot")),
      _displacement_jump_global_dot(
          declareProperty<RealVectorValue>("displacement_jump_global_dot")),
      _disp_dot(_ndisp), _disp_neighbor_dot(_ndisp)

{
  // initializing the displacement vectors
  for (unsigned int i = 0; i < _ndisp; ++i) {
    _disp_dot[i] = &coupledDot("displacements", i);
    _disp_neighbor_dot[i] = &coupledNeighborValueDot("displacements", i);
  }
}

void ViscousSlidingCZM::computeQpProperties() {

  RealTensorValue RotationGlobalToLocal =
      RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));

  // computing the displacement jump rate
  for (unsigned int i = 0; i < _ndisp; i++)
    _displacement_jump_global_dot[_qp](i) =
        (*_disp_neighbor_dot[i])[_qp] - (*_disp_dot[i])[_qp];
  for (unsigned int i = _ndisp; i < 3; i++)
    _displacement_jump_global_dot[_qp](i) = 0;

  _displacement_jump_dot[_qp] =
      RotationGlobalToLocal * _displacement_jump_global_dot[_qp];

  PureElasticCZM::computeQpProperties();
}

void ViscousSlidingCZM::initQpStatefulProperties() {
  for (unsigned int i = 0; i < 3; i++)
    _traction[_qp](i) = 0;
}

void ViscousSlidingCZM::ComputeShearTraction(RealVectorValue &traction) {
  Real C = ComputeShearStiffness();
  for (unsigned int i = 1; i < 3; i++) {
    Real eta = _shear_viscosity;
    traction(i) =
        eta * _displacement_jump_dot[_qp](i) +
        std::exp(-_dt * C / eta) *
            (_traction_old[_qp](i) - eta * _displacement_jump_dot[_qp](i));
  }
}

void ViscousSlidingCZM::ComputeShearTractionDerivatives(
    RankTwoTensor &traction_derivatives) {
  Real C = ComputeShearStiffness();
  RankTwoTensor dC_dui = ComputeShearStiffnessDerivatives();
  for (unsigned int i = 1; i < 3; i++) {
    Real eta = _shear_viscosity;
    Real dTsi_dui = eta * (1. - std::exp(-_dt * C / eta)) / _dt;
    Real dTsi_dC =
        -_dt / eta * std::exp(-_dt * C / eta) *
        (_traction_old[_qp](i) - eta * _displacement_jump_dot[_qp](i));

    traction_derivatives(i, i) = dTsi_dui + dTsi_dC * dC_dui(i, i);
  }
}
