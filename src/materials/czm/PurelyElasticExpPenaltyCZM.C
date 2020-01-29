//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceValueTools.h"
#include "PurelyElasticExpPenaltyCZM.h"

registerMooseObject("DeerApp", PurelyElasticExpPenaltyCZM);

template <> InputParameters validParams<PurelyElasticExpPenaltyCZM>() {
  InputParameters params = validParams<CZMMaterialBase>();
  params.addClassDescription(
      "time dependent interface damage model with linear interpolation");
  params.addRequiredParam<Real>("K_n", "stiffness in the normal direction");
  params.addRequiredParam<Real>("K_t", "stiffness in the shear direction");
  params.addParam<Real>("copenetration_penalty_linear", 100,
                        "the copenetration penalty for bilinear method");
  params.addParam<Real>("exponential_penalty_multiplier", 1e4,
                        "the copenetration penalty");
  params.addParam<Real>("exponential_shift", 1e-3,
                        "the dispalcement jump at which penalty kicks in");
  params.addParam<Real>("exponential_characteristicopening_length", 1e-6,
                        "the dispalcement jump at which penalty kicks in");
  return params;
}

PurelyElasticExpPenaltyCZM::PurelyElasticExpPenaltyCZM(
    const InputParameters &parameters)
    : CZMMaterialBase(parameters),
      _stiffness({getParam<Real>("K_n"), getParam<Real>("K_t"),
                  getParam<Real>("K_t")}),
      _copenetration_penalty(getParam<Real>("copenetration_penalty_linear")),
      _exp_penalty(getParam<Real>("exponential_penalty_multiplier")),
      _exp_shift(getParam<Real>("exponential_shift")),
      _exp_length(getParam<Real>("exponential_characteristicopening_length"))

{}

RealVectorValue PurelyElasticExpPenaltyCZM::computeTraction() {
  // The convention for ordering the traction is N, T, S, where N is the normal
  // direction, and T and S are two arbitrary tangential directions.
  RealVectorValue traction_local;

  for (unsigned int i = 0; i < 3; i++)
    traction_local(i) = _stiffness[i] * _displacement_jump[_qp](i);

  traction_local(0) += TractionPenaltyExponential();
  return traction_local;
}

RankTwoTensor PurelyElasticExpPenaltyCZM::computeTractionDerivatives() {
  RankTwoTensor traction_jump_derivatives_local;

  for (unsigned int i = 0; i < 3; i++)
    traction_jump_derivatives_local(i, i) = _stiffness[i];

  traction_jump_derivatives_local(0, 0) +=
      TractionPenaltyExponentialDerivative();
  return traction_jump_derivatives_local;
}

Real PurelyElasticExpPenaltyCZM::TractionPenaltyExponential() const {

  return -_copenetration_penalty *
             std::exp(-(_displacement_jump[_qp](0) / _exp_length + _exp_shift) *
                      _exp_penalty) +
         _copenetration_penalty * std::exp(_exp_shift) * _exp_penalty;
}

Real PurelyElasticExpPenaltyCZM::TractionPenaltyExponentialDerivative() const {
  return _copenetration_penalty *
         std::exp(-(_displacement_jump[_qp](0) / _exp_length + _exp_shift) *
                  _exp_penalty) *
         _exp_penalty / _exp_length;
}
