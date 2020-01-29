//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceValueTools.h"
#include "TimeDependentDamage.h"

registerMooseObject("DeerApp", TimeDependentDamage);

template <> InputParameters validParams<TimeDependentDamage>() {
  InputParameters params = validParams<CZMMaterialBase>();
  params.addClassDescription(
      "time dependent interface damage model with linear interpolation");
  params.addRequiredParam<Real>("K_n", "stiffness in the normal direction");
  params.addRequiredParam<Real>("K_t", "stiffness in the shear direction");
  params.addParam<Real>("copenetration_penalty", 1e2,
                        "the copenetration penalty");
  params.addParam<Real>("stiffness_reduction_factor", 100,
                        "the stiffness reduction factor used to "
                        "determine the minimum rallowed stffness");
  params.addParam<Real>("max_damage", 0.95,
                        "the maximum damage before element becomes inactive");

  params.addParam<Real>("fluid_pressure", 0,
                        "the fluid pressure to apply after elemnt failure");

  params.addRequiredParam<MaterialPropertyName>(
      "effective_stress_mp_name",
      "the effective stress material property name");
  params.addParam<Real>(
      "min_allowed_residual_life", 0.,
      "the minimum residaul life before traction decay kicks in");
  // params.addParam<MooseEnum>("interface_value_type",
  //                            InterfaceValueTools::InterfaceAverageOptions(),
  //                            "effective stress averaging type");
  params.addParam<std::vector<Real>>(
      "x", "The abscissa values for interpolating the damage (e.g stress)");
  params.addParam<std::vector<Real>>(
      "y", "The ordinate values for interpolating the damage (e.g. life)");
  params.addParam<Real>("residual_life_scaling_factor", 1.,
                        "How quickly the traction should be forced to zero");
  return params;
}

TimeDependentDamage::TimeDependentDamage(const InputParameters &parameters)
    : CZMMaterialBase(parameters),
      _stiffness({getParam<Real>("K_n"), getParam<Real>("K_t"),
                  getParam<Real>("K_t")}),
      _minimum_stiffnes(
          {_stiffness[0] / getParam<Real>("stiffness_reduction_factor"),
           _stiffness[1] / getParam<Real>("stiffness_reduction_factor"),
           _stiffness[2] / getParam<Real>("stiffness_reduction_factor")}),
      _copenetration_penalty(getParam<Real>("copenetration_penalty")),
      _max_damage(getParam<Real>("max_damage")),
      _fluid_pressure(getParam<Real>("fluid_pressure")),
      _damage(declareProperty<Real>("interface_damage")),
      _damage_old(getMaterialPropertyOld<Real>("interface_damage")),
      _element_failed(declareProperty<Real>("element_failed")),
      _element_failed_old(getMaterialPropertyOld<Real>("element_failed")),
      _effective_stress_old(getMaterialPropertyOldByName<Real>(
          getParam<MaterialPropertyName>("effective_stress_mp_name"))),
      _residual_life(declareProperty<Real>("residual_life")),
      _time_fail(declareProperty<Real>("failure_time")),
      _du_fail(declareProperty<RealVectorValue>("jump_at_failure")),
      _T_fail(declareProperty<RealVectorValue>("traction_at_failure")),
      _K_fail(declareProperty<RealVectorValue>("stiffness_at_failure")),
      _residual_life_old(getMaterialPropertyOld<Real>("residual_life")),
      _time_fail_old(getMaterialPropertyOld<Real>("failure_time")),
      _du_fail_old(getMaterialPropertyOld<RealVectorValue>("jump_at_failure")),
      _T_fail_old(
          getMaterialPropertyOld<RealVectorValue>("traction_at_failure")),
      _K_fail_old(
          getMaterialPropertyOld<RealVectorValue>("stiffness_at_failure")),
      _residual_life_scaling_factor(
          getParam<Real>("residual_life_scaling_factor")),
      _min_allowed_residual_life(getParam<Real>("min_allowed_residual_life"))

{
  // check stifness reduction factor
  if (getParam<Real>("stiffness_reduction_factor") <= 1.)
    mooseError("TimeDependentDamage: stiffness_reduction_factor must be grater "
               "than 1");

  for (unsigned int i = 0; i < 3; i++)
    if (_stiffness[i] * (1. - _max_damage) < _minimum_stiffnes[i])
      mooseError("TimeDependentDamage: stiffness_reduction_factor and "
                 "max_damage are conflicting as stiffness*(1-max_damage)< "
                 "stiffness/stiffness_reduction_factor");

  std::vector<Real> x;
  std::vector<Real> y;

  if ((parameters.isParamValid("x")) || (parameters.isParamValid("y"))) {
    if (!((parameters.isParamValid("x")) && (parameters.isParamValid("y"))))
      mooseError(
          "In PiecewiseLinearInterpolationMaterial ", _name,
          ": Both 'x' and 'y' must be specified if either one is specified.");

    if (parameters.isParamValid("xy_data"))
      mooseError("In PiecewiseLinearInterpolationMaterial ", _name,
                 ": Cannot specify 'x', 'y', and 'xy_data' together.");

    x = getParam<std::vector<Real>>("x");
    y = getParam<std::vector<Real>>("y");
  }

  try {
    _linear_interp = libmesh_make_unique<LinearInterpolation>(x, y);
  } catch (std::domain_error &e) {
    mooseError("In TimeDependentDamage ", _name, ": ", e.what());
  }
}

RealVectorValue TimeDependentDamage::computeTraction() {
  // The convention for ordering the traction is N, T, S, where N is the normal
  // direction, and T and S are two arbitrary tangential directions.
  RealVectorValue traction_local;

  // update damage and check if element has failed
  updateDamage();
  if (_element_failed_old[_qp] == 0) {
    for (unsigned int i = 0; i < 3; i++)
      traction_local(i) =
          _stiffness[i] * (1. - _damage[_qp]) * _displacement_jump[_qp](i);

    computeDecayRelatedVariables(traction_local);
  } else
    traction_local = computeExpTractionDecay();

  if (_displacement_jump[_qp](0) < 0)
    traction_local(0) *= _copenetration_penalty;

  if (_element_failed_old[_qp] > 0)
    traction_local(0) -= _fluid_pressure;

  _element_failed[_qp] = 0;
  if ((_element_failed_old[_qp] == 1) ||
      ((_element_failed_old[_qp] == 0) &&
       ((_damage[_qp] >= _max_damage) ||
        _residual_life[_qp] < _min_allowed_residual_life)))
    _element_failed[_qp] = 1;

  return traction_local;
}

RankTwoTensor TimeDependentDamage::computeTractionDerivatives() {
  RankTwoTensor traction_jump_derivatives_local;

  if (_element_failed_old[_qp] == 0) {
    for (unsigned int i = 0; i < 3; i++)
      traction_jump_derivatives_local(i, i) =
          _stiffness[i] * (1. - _damage[_qp]);

  } else
    traction_jump_derivatives_local = computeExpTractionDerivativeDecay();

  if (_displacement_jump[_qp](0) < 0)
    traction_jump_derivatives_local(0, 0) *= _copenetration_penalty;
  return traction_jump_derivatives_local;
}

void TimeDependentDamage::updateDamage() {

  const Real max_life = _linear_interp->sample(_effective_stress_old[_qp]);
  Real delta_damage = _dt / max_life;

  _damage[_qp] = std::min(_damage_old[_qp] + delta_damage, _max_damage);
}

void TimeDependentDamage::initQpStatefulProperties() {
  _damage[_qp] = 0;
  _element_failed[_qp] = 0;
  _time_fail[_qp] = 0;
  _residual_life[_qp] = _linear_interp->sample(0.);
  for (unsigned int i = 0; i < 3; i++) {
    _du_fail[_qp](i) = 0;
    _T_fail[_qp](i) = 0;
    _K_fail[_qp](i) = _stiffness[i];
  }
}

void TimeDependentDamage::computeDecayRelatedVariables(
    const RealVectorValue &traction_local) {
  if (_element_failed_old[_qp] == 0) { // update values
    const Real max_life = _linear_interp->sample(_effective_stress_old[_qp]);
    _residual_life[_qp] =
        std::max((1. - _damage[_qp]) * max_life, _min_allowed_residual_life);
    _time_fail[_qp] = _t;
    _T_fail[_qp] = traction_local;
    _du_fail[_qp] = _displacement_jump[_qp];
    for (unsigned int i = 0; i < 3; i++)
      _K_fail[_qp](i) = _stiffness[i] * (1. - _damage[_qp]);
  } else {
    _residual_life[_qp] = _residual_life_old[_qp];
    _time_fail[_qp] = _time_fail_old[_qp];
    _T_fail[_qp] = _T_fail_old[_qp];
    _du_fail[_qp] = _du_fail_old[_qp];
    _K_fail[_qp] = _K_fail_old[_qp];
  }
}

RealVectorValue TimeDependentDamage::computeExpTractionDecay() {

  RealVectorValue traction_local;
  computeDecayRelatedVariables(traction_local);

  Real damage_exp = std::exp(-5. * _residual_life_scaling_factor *
                             (_t - _time_fail[_qp]) / _residual_life[_qp]);

  for (unsigned int i = 0; i < 3; i++) {
    Real c_exp = std::max(_K_fail[_qp](i) * damage_exp, _minimum_stiffnes[i]);

    traction_local(i) =
        (_displacement_jump[_qp](i) - _du_fail[_qp](i)) * c_exp +
        _T_fail[_qp](i) * damage_exp;
  }
  return traction_local;
}

RankTwoTensor TimeDependentDamage::computeExpTractionDerivativeDecay() {
  RankTwoTensor traction_jump_derivatives_local;
  Real damage_exp = std::exp(-5. * _residual_life_scaling_factor *
                             (_t - _time_fail[_qp]) / _residual_life[_qp]);
  for (unsigned int i = 0; i < 3; i++) {

    Real c_exp = std::max(_K_fail[_qp](i) * damage_exp, _minimum_stiffnes[i]);
    traction_jump_derivatives_local(i, i) = c_exp;
    // std::cout << "traction_jump_derivatives_local(" << i << "," << i
    //           << ")= " << traction_jump_derivatives_local(i, i) << std::endl;
  }
  return traction_jump_derivatives_local;
}
