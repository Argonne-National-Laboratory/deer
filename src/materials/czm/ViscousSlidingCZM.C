//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViscousSlidingCZM.h"

registerMooseObject("DeerApp", ViscousSlidingCZM);

InputParameters
ViscousSlidingCZM::validParams()
{
  InputParameters params = PureElasticCZM::validParams();
  params.addClassDescription("Cohesive model with linear leastic opening and viscous sliding");
  params.addRequiredParam<Real>("shear_viscosity", "Interface shear viscosity");

  return params;
}

ViscousSlidingCZM::ViscousSlidingCZM(const InputParameters & parameters)
  : PureElasticCZM(parameters),
    _shear_viscosity(getParam<Real>("shear_viscosity")),
    _interface_traction_old(
        getMaterialPropertyOldByName<RealVectorValue>(_base_name + "interface_traction")),
    _interface_displacement_jump_old(
        getMaterialPropertyOldByName<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_displacement_jump_dot(
        declarePropertyByName<RealVectorValue>(_base_name + "interface_displacement_jump_dot"))
{
}

void
ViscousSlidingCZM::computeQpProperties()
{

  _interface_displacement_jump_dot[_qp] =
      (_interface_displacement_jump[_qp] - _interface_displacement_jump_old[_qp]) / _dt;

  PureElasticCZM::computeQpProperties();
}

void
ViscousSlidingCZM::initQpStatefulProperties()
{
  PureElasticCZM::initQpStatefulProperties();
  _interface_traction[_qp] = 0;
}

void
ViscousSlidingCZM::ComputeShearTraction(RealVectorValue & traction)
{
  Real C = ComputeShearStiffness();
  Real eta = ComputeShearViscosity();
  for (unsigned int i = 1; i < 3; i++)
  {
    traction(i) = eta * _interface_displacement_jump_dot[_qp](i) +
                  std::exp(-_dt * C / eta) * (_interface_traction_old[_qp](i) -
                                              eta * _interface_displacement_jump_dot[_qp](i));
  }
}

void
ViscousSlidingCZM::ComputeShearTractionDerivatives(RankTwoTensor & traction_derivatives)
{
  Real C = ComputeShearStiffness();
  RankTwoTensor dC_dui = ComputeShearStiffnessDerivatives();
  Real eta = ComputeShearViscosity();
  RankTwoTensor deta_dui = ComputeShearViscosityDerivatives();
  for (unsigned int i = 1; i < 3; i++)
  {
    Real dTsi_dui = eta * (1. - std::exp(-_dt * C / eta)) / _dt;
    Real dTsi_dC =
        -_dt / eta * std::exp(-_dt * C / eta) *
        (_interface_traction_old[_qp](i) - eta * _interface_displacement_jump_dot[_qp](i));
    Real dTsi_deta =
        _interface_displacement_jump_dot[_qp](i) -
        _dt * C / (eta * eta) * std::exp(-_dt * C / eta) *
            (_interface_traction_old[_qp](i) - eta * _interface_displacement_jump_dot[_qp](i)) -
        std::exp(-_dt * C / eta) * _interface_displacement_jump_dot[_qp](i);

    traction_derivatives(i, i) = dTsi_dui + dTsi_dC * dC_dui(i, i) + dTsi_deta * deta_dui(i, i);
  }
}

Real
ViscousSlidingCZM::ComputeShearViscosity()
{
  return _shear_viscosity;
}

RankTwoTensor
ViscousSlidingCZM::ComputeShearViscosityDerivatives()
{
  // in this case derivatives are all zero so we just retrun a zero tensor;
  return RankTwoTensor::initNone;
}
