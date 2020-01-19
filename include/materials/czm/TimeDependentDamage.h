//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CZMMaterialBase.h"
#include "LinearInterpolation.h"
class TimeDependentDamage;
template <> InputParameters validParams<TimeDependentDamage>();
/**
 * Implementation of the non-stateful exponential traction separation law
 * proposed by Salehani, Mohsen Khajeh and Irani, Nilgoon 2018
 **/
class TimeDependentDamage : public CZMMaterialBase {
public:
  TimeDependentDamage(const InputParameters &parameters);

protected:
  virtual void initQpStatefulProperties() override;

  virtual void updateDamage();

  virtual RealVectorValue computeTraction() override;

  virtual RankTwoTensor computeTractionDerivatives() override;

  /// the displacement jump associated to the maximum traction
  const std::vector<Real> _stiffness;
  const std::vector<Real> _minimum_stiffnes;
  const Real _copenetration_penalty;
  const Real _max_damage;
  const Real _fluid_pressure;
  MaterialProperty<Real> &_damage;
  const MaterialProperty<Real> &_damage_old;
  MaterialProperty<Real> &_element_failed;
  const MaterialProperty<Real> &_element_failed_old;
  const MaterialProperty<Real> &_effective_stress_old;

  void computeDecayRelatedVariables(const RealVectorValue &traction_local);
  RealVectorValue computeExpTractionDecay();
  RankTwoTensor computeExpTractionDerivativeDecay();
  MaterialProperty<Real> &_residual_life;
  MaterialProperty<Real> &_time_fail;
  MaterialProperty<RealVectorValue> &_du_fail;
  MaterialProperty<RealVectorValue> &_T_fail;
  MaterialProperty<RealVectorValue> &_K_fail;
  const MaterialProperty<Real> &_residual_life_old;
  const MaterialProperty<Real> &_time_fail_old;
  const MaterialProperty<RealVectorValue> &_du_fail_old;
  const MaterialProperty<RealVectorValue> &_T_fail_old;
  const MaterialProperty<RealVectorValue> &_K_fail_old;
  const Real _residual_life_scaling_factor;

  /// LinearInterpolation object
  std::unique_ptr<LinearInterpolation> _linear_interp;
};
