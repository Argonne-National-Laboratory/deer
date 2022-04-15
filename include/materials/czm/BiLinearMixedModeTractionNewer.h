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
class BiLinearMixedModeTractionNewer : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  BiLinearMixedModeTractionNewer(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeInterfaceTractionAndDerivatives() override;

  /// Calculate the traction and derivative
  virtual std::tuple<RealVectorValue, RankTwoTensor> updateState(const
                                                                 RealVectorValue
                                                                 & delta);
  
  /// calculate the mode mixity
  virtual std::tuple<Real,RealVectorValue> modeMixity(const RealVectorValue & delta) const;

  /// calculate the damage
  virtual std::tuple<Real,Real,Real,Real> damage(const Real & delta, const Real & delta_init, 
                                                 const Real & delta_final) const;

  /// calculate the effective displacement
  std::tuple<Real, RealVectorValue> deltaEffective(const RealVectorValue & delta) const;

  /// calculate the onset and final displacements
  virtual std::tuple<Real,Real,Real,Real> deltaThreshold(const RealVectorValue & delta, 
                                                          const Real & beta) const;

  /// Old displacements (for lagging)
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;

  /// penalty elastic stiffness
  const Real _K;

  ///@{
  /// damage variable
  MaterialProperty<Real> & _d;
  const MaterialProperty<Real> & _d_old;
  ///@}
  
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

  /// mixed mode propagation criterion
  enum class MixedModeCriterion
  {
    POWER_LAW,
    BK
  } _criterion;

  /// Lag the calculation of the mode-mixity
  bool _lag_beta;

  /// Lag the calculation of the damage
  bool _lag_damage;
};
