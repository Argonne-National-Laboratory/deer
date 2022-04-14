//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CZMComputeLocalTractionTotalBase.h"

/**
 * Implementation of the mixed mode bilinear traction separation law
 * described in Mixed-Mode Decohesion Finite Elements for the Simulation of Delamination in
 *Composite Materials, Pedro P. Camanho and Carlos G. Davila, NASA/TM-2002-211737
 **/
class BiLinearMixedModeTractionNew : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  BiLinearMixedModeTractionNew(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeInterfaceTractionAndDerivatives() override;

  virtual Real modeMixity(const RealVectorValue & delta);

  virtual Real damage(const Real & delta, const Real & delta_init, const Real & delta_final);

  /// method computing the total traction
  virtual RealVectorValue computeTraction();

  /// method computing the total traction derivatives w.r.t. the interface displacement jump
  virtual RankTwoTensor computeTractionDerivatives();

  /// penalty elastic stiffness
  const Real _K;

  ///@{
  /// damage variable
  MaterialProperty<Real> & _d;
  const MaterialProperty<Real> & _d_old;
  ///@}

  /// old interface displacement jump value
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;

  ///@{
  /// critical Mode I and II fracture toughness
  const Real _GI_C;
  const Real _GII_C;
  ///@}

  ///@{
  /// onset normal seperation and shear seperation
  const Real _N;
  const Real _S;
  ///@}

  /// The B-K power law parameter
  const Real _eta;

  /// mode_mixity_ratio
  MaterialProperty<Real> & _beta;

  /// mixed mode propagation criterion
  enum class MixedModeCriterion
  {
    POWER_LAW,
    BK
  } _criterion;
};
