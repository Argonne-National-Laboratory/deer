//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PureElasticCZM.h"

registerMooseObject("DeerApp", PureElasticCZM);

InputParameters
PureElasticCZM::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Cohesive model with linear elastic opening and shearing");
  params.addRequiredParam<Real>("E", "Interface normal elastic modulus");
  params.addRequiredParam<Real>("G", "Interface shear elastic modulus");
  params.addRequiredParam<Real>("interface_thickness", "Initial interface thickness");
  params.addParam<Real>("penetration_penalty",
                        1e2,
                        "The penetration penalty applied after (duN < -interface_thickness)");
  return params;
}

PureElasticCZM::PureElasticCZM(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _E(getParam<Real>("E")),
    _G(getParam<Real>("G")),
    _interface_thickness(getParam<Real>("interface_thickness")),
    _penetration_penalty(getParam<Real>("penetration_penalty"))
{
}

void
PureElasticCZM::computeInterfaceTractionAndDerivatives()
{
  _interface_traction[_qp] = computeTraction();
  _dinterface_traction_djump[_qp] = computeTractionDerivatives();
}

RealVectorValue
PureElasticCZM::computeTraction()
{

  RealVectorValue traction;

  ComputeNormalTraction(traction);
  ComputeShearTraction(traction);

  return traction;
}

RankTwoTensor
PureElasticCZM::computeTractionDerivatives()
{
  RankTwoTensor traction_jump_derivatives;

  ComputeNormalTractionDerivatives(traction_jump_derivatives);
  ComputeShearTractionDerivatives(traction_jump_derivatives);

  return traction_jump_derivatives;
}

void
PureElasticCZM::ComputeShearTraction(RealVectorValue & traction)
{
  // note that according to the CZM implementaion in MOOSE we always work with
  // 3D objects for disaplcement jump, traction and derivatives.
  // The rotation takes care of the rest
  Real C = ComputeShearStiffness();
  for (unsigned int i = 1; i < 3; i++)
    traction(i) = C * _interface_displacement_jump[_qp](i);
}

void
PureElasticCZM::ComputeShearTractionDerivatives(RankTwoTensor & traction_derivatives)
{
  Real C = ComputeShearStiffness();
  RankTwoTensor dstiffness_du = ComputeShearStiffnessDerivatives();

  for (unsigned int i = 1; i < 3; i++)
  {
    traction_derivatives(i, i) = C;
    for (unsigned int j = 0; i < 3; i++)
      traction_derivatives(i, j) += _interface_displacement_jump[_qp](i) * dstiffness_du(i, j);
  }
}

void
PureElasticCZM::ComputeNormalTraction(RealVectorValue & traction)
{
  Real C = ComputeNormalStiffness();
  traction(0) = C * _interface_displacement_jump[_qp](0);

  if (_interface_displacement_jump[_qp](0) < -_interface_thickness)
    traction(0) *= _penetration_penalty;
}

void
PureElasticCZM::ComputeNormalTractionDerivatives(RankTwoTensor & traction_derivatives)
{

  Real C = ComputeNormalStiffness();
  RealVectorValue dstiffness_du = ComputeNormalStiffnessDerivatives();

  traction_derivatives(0, 0) = C;
  for (unsigned int i = 0; i < 3; i++)
    traction_derivatives(0, i) += _interface_displacement_jump[_qp](0) * dstiffness_du(i);

  if (_interface_displacement_jump[_qp](0) < -_interface_thickness)
    traction_derivatives(0, 0) *= _penetration_penalty;
}

Real
PureElasticCZM::ComputeNormalStiffness()
{
  return _E / _interface_thickness;
}

RealVectorValue
PureElasticCZM::ComputeNormalStiffnessDerivatives()
{
  RealVectorValue dstiffness_du;
  return dstiffness_du;
}

Real
PureElasticCZM::ComputeShearStiffness()
{
  return _G / _interface_thickness;
}

RankTwoTensor
PureElasticCZM::ComputeShearStiffnessDerivatives()
{
  RankTwoTensor dstiffness_dui;
  return dstiffness_dui;
}
