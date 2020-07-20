//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CZMMaterialBasePK1.h"
#include "ShamNeedlemanEquation.h"
/**
 * Implementation of teh grain boundary cavitation model proposed by Sham
 * and Needlemena 1983 and expanded by Van der Gieseen 1995 **/
class GBCavitation : public CZMMaterialBasePK1 {
public:
  static InputParameters validParams();
  GBCavitation(const InputParameters &parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeTractionIncrementAndDerivatives() override;
  void computeAverageBulkPorperties();
  void initNLSystemParamter(std::vector<std::string> &pname, vecD &pvalue,
                            std::vector<std::string> &rate_pname,
                            vecD &rate_pvalue);

  /// method to kill the traction
  void tractionDecay();

  const bool _use_old_bulk_property;
  /// reference to bulk properties
  ///@{
  const MaterialProperty<RankTwoTensor> &_stress_master;
  const MaterialProperty<RankTwoTensor> &_stress_slave;
  const MaterialProperty<RankTwoTensor> &_inelastic_strain_master;
  const MaterialProperty<RankTwoTensor> &_inelastic_strain_slave;
  const MaterialProperty<RankTwoTensor> &_inelastic_strain_master_old;
  const MaterialProperty<RankTwoTensor> &_inelastic_strain_slave_old;
  ///@}

  /// the computed equivalent stress and strain values on the interface
  ///@{
  MaterialProperty<Real> &_stress_vm;
  const MaterialProperty<Real> &_stress_vm_old;
  MaterialProperty<Real> &_stress_H;
  const MaterialProperty<Real> &_stress_H_old;
  MaterialProperty<Real> &_strain_rate_eq;
  const MaterialProperty<Real> &_strain_rate_eq_old;
  MaterialProperty<Real> &_strain_eq;
  const MaterialProperty<Real> &_strain_eq_old;
  ///@}

  /// sham needleman state variables
  ///@{
  MaterialProperty<Real> &_a;
  const MaterialProperty<Real> &_a_old;
  MaterialProperty<Real> &_b;
  const MaterialProperty<Real> &_b_old;
  MaterialProperty<int> &_nucleation_is_active;
  const MaterialProperty<int> &_nucleation_is_active_old;
  MaterialProperty<Real> &_D;
  MaterialProperty<Real> &_D_rate;
  ///@}

  /// failure state variables
  ///@{
  MaterialProperty<int> &_element_failed;
  const MaterialProperty<int> &_element_failed_old;
  MaterialProperty<Real> &_time_at_failure;
  const MaterialProperty<Real> &_time_at_failure_old;
  MaterialProperty<RealVectorValue> &_traction_at_failure;
  const MaterialProperty<RealVectorValue> &_traction_at_failure_old;
  MaterialProperty<RealVectorValue> &_jump_at_failure;
  const MaterialProperty<RealVectorValue> &_jump_at_failure_old;
  MaterialProperty<Real> &_residual_life;
  const MaterialProperty<Real> &_residual_life_old;

  ///@}

  /// sham needleman eqauations paramters
  ///@{
  const Real _a0;
  const Real _b0;
  const Real _NI;
  const Real _FN_NI;
  const Real _FN;
  const Real _Nmax_NI;
  const Real _b_sat;
  const Real _S0;
  const Real _beta;
  const Real _psi_degree;
  const Real _h;
  const Real _E_GB;
  const Real _G_GB;
  const Real _D_GB;
  const Real _eta_sliding;
  const Real _thickness;
  const Real _n;
  const Real _theta;
  const Real _E_penalty_minus_thickenss_over_2;
  const Real _E_penalty_minus_thickenss;
  ///@}

  /// switch for activing physics
  ///@{
  const bool _nucleation_on;
  const bool _growth_on;
  const bool _use_triaxial_growth;
  ///@}

  /// failure constants
  ///@{
  const Real _D_failure;
  const Real _minimum_allowed_residual_life;
  const Real _maximum_allowed_opening_traction;
  const Real _minimum_allowed_stiffness;
  ///@}

  /// non linear solver paramters
  ///@{
  const unsigned int _max_time_cut;
  const unsigned int _max_nonlinear_iter;
  const Real _nl_residual_abs_tol;
  const bool _force_substep;
  ///@}
};
