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

class PurelyElasticExpPenaltyCZM;
template <> InputParameters validParams<PurelyElasticExpPenaltyCZM>();
/**
 * Implementation of the non-stateful exponential traction separation law
 * proposed by Salehani, Mohsen Khajeh and Irani, Nilgoon 2018
 **/
class PurelyElasticExpPenaltyCZM : public CZMMaterialBase {
public:
  PurelyElasticExpPenaltyCZM(const InputParameters &parameters);

protected:
  virtual RealVectorValue computeTraction() override;

  virtual RankTwoTensor computeTractionDerivatives() override;

  Real TractionPenaltyExponential() const;
  Real TractionPenaltyExponentialDerivative() const;

  /// the displacement jump associated to the maximum traction
  const std::vector<Real> _stiffness;
  const Real _copenetration_penalty;
  const Real _exp_penalty;
  const Real _exp_shift;
  const Real _exp_length;
};
