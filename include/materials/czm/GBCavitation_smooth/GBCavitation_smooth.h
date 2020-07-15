//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "petscmat.h"
#include "petscsnes.h"

#include "CZMMaterialBasePK1.h"

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

#include "GBCavitationNLSystem.h"
#include "scaleFactors.h"

#include "GBCavitationBoundaryPropertyUO.h"

/* DATA STRUCTURE TO PASS QUANTITIES TO PETSC*/
typedef struct {
  PetscReal dt;
  GBCavitationNLSystem *GBNLsysstem_pt;
  const nlFunBase::io_maps_type *scale_factor_pt;
  nlFunBase::io_maps_type *my_xold;
  nlFunBase::io_maps_type *my_params;
  const int *n_eq_pt;
} ApplicationCtx;

/**
 * This material automatically declares as material properties whatever is
 * passed to it through the parameters 'prop_names' and uses the values from
 * 'prop_values' as the values for those properties.
 *
 * This is not meant to be used in a production capacity... and instead is meant
 * to be used during development phases for ultimate flexibility.
 */
class GBCavitation_smooth : public CZMMaterialBasePK1 {
public:
  static InputParameters validParams();
  GBCavitation_smooth(const InputParameters &parameters);
  ~GBCavitation_smooth();

protected:
  const GBCavitationBoundaryPropertyUO *_GBCavitationBoundaryPropertyUO;
  MaterialProperty<RealVectorValue> &_displacement_jump_dot;
  const MaterialProperty<RealVectorValue> &_displacement_jump_dot_old;

  /* INTERFACE PARAMETERS */
  // const Real _a0;                  /*initial cavity half radius*/
  // const Real _b0;                  /*initial cavity half spacing*/
  // const Real _b_saturation;        /*saturation cavity half spacing*/
  // const Real _NI;                  /*number of intial cavity*/
  // const Real _FN;                  /*cavitration rate*/
  // const Real _S_thr;               /*thresould value for cavitation to
  // occur*/ const Real _interface_thickness; /*interface Young modulus*/ const
  // Real _E_interface;         /*interface Young modulus*/
  const Real _E_penalty; /*incase of copenetration*/
  // const Real _G_interface;
  // /*interface shear modulus*/
  const Real _beta_exponent; /*Traction exponent*/
  // const Real _D_gb;          /*gran boundary diffusion
  // coefficient*/
  const Real _n_exponent; /*power law creep exponent*/
  const Real _alpha_n;    /*3/(2*_n_exponent))*/
  // const Real _psi_angle;   /*equilibrium cavity tip half-angle [degree]*/
  // const Real _h;           /*function of psi angle*/
  // const Real _sigma_0;     /*traction normalization parameter*/
  // const Real _eta_sliding; /*interface sliding viscosity*/

  void getInitPropertyValuesFromParams(Real &FN_NI, Real &Nmax_NI, Real &a0,
                                       Real &b0, Real &psi, Real &D_gb,
                                       Real &E_interface, Real &G_interface,
                                       Real &eta_sliding, Real &sigma_0,
                                       Real &S_thr) const;

  void getInitPropertyValuesFromUO(Real &FN_NI, Real &Nmax_NI, Real &a0,
                                   Real &b0, Real &psi, Real &D_gb,
                                   Real &E_interface, Real &G_interface,
                                   Real &eta_sliding, Real &sigma_0,
                                   Real &S_thr) const;

  void InitGBCavitationParamsAndProperties();

  void computeAverageBulkPorperties();

  MaterialProperty<Real> &_a0;
  const MaterialProperty<Real> &_a0_old;

  MaterialProperty<Real> &_b0;
  const MaterialProperty<Real> &_b0_old;

  MaterialProperty<Real> &_NI;
  const MaterialProperty<Real> &_NI_old;

  MaterialProperty<Real> &_FN;
  const MaterialProperty<Real> &_FN_old;

  MaterialProperty<Real> &_D_gb;
  const MaterialProperty<Real> &_D_gb_old;

  MaterialProperty<Real> &_b_saturation;
  const MaterialProperty<Real> &_b_saturation_old;

  MaterialProperty<Real> &_E_interface;
  const MaterialProperty<Real> &_E_interface_old;

  MaterialProperty<Real> &_G_interface;
  const MaterialProperty<Real> &_G_interface_old;

  MaterialProperty<Real> &_eta_sliding;
  const MaterialProperty<Real> &_eta_sliding_old;

  MaterialProperty<Real> &_interface_thickness;
  const MaterialProperty<Real> &_interface_thickness_old;

  MaterialProperty<Real> &_sigma_0;
  const MaterialProperty<Real> &_sigma_0_old;

  MaterialProperty<Real> &_S_thr;
  const MaterialProperty<Real> &_S_thr_old;

  MaterialProperty<Real> &_h;
  const MaterialProperty<Real> &_h_old;

  // GBCavitationNLSystem _GBNLsystem;

  /* helper fucntion to translate petsc variable to interface inputs*/
  scaleFactors _residualScaleFactors;
  nlFunBase::io_maps_type _residual_scale_factors;

  // virtual RealVectorValue computeTraction() override {
  //   RealVectorValue a;
  //   return a;
  // };
  // virtual RankTwoTensor computeTractionDerivatives() override {
  //   RankTwoTensor a;
  //   return a;
  // };

  virtual void computeTractionIncrementAndDerivatives() override;
  void initQpStatefulProperties() override;

  // required additional material properties
  MaterialProperty<Real> &_a;
  const MaterialProperty<Real> &_a_old;
  MaterialProperty<Real> &_b;
  const MaterialProperty<Real> &_b_old;
  MaterialProperty<Real> &_D;
  const MaterialProperty<Real> &_D_old;
  MaterialProperty<Real> &_D_dot;
  const MaterialProperty<Real> &_D_dot_old;
  MaterialProperty<Real> &_residual_life;
  const MaterialProperty<Real> &_residual_life_old;
  MaterialProperty<Real> &_VL1_dot;
  MaterialProperty<Real> &_VL2_dot;
  MaterialProperty<Real> &_VH1_dot;
  MaterialProperty<Real> &_VH2_dot;
  MaterialProperty<Real> &_Vdot;
  MaterialProperty<Real> &_VLdot;
  const MaterialProperty<Real> &_VLdot_old;
  MaterialProperty<Real> &_VHdot;
  const MaterialProperty<Real> &_VHdot_old;
  MaterialProperty<bool> &_nucleation_above_threshold;
  const MaterialProperty<bool> &_nucleation_above_threshold_old;
  /* BULK MATERIAL PROPERTIES REQUIRED FOR INTERFACE CALCUALTION */
  const bool _use_old_avg_prop;
  /// reference to bulk properties
  ///@{
  const MaterialProperty<RankTwoTensor> &_stress_master;
  const MaterialProperty<RankTwoTensor> &_stress_slave;
  const MaterialProperty<RankTwoTensor> &_inelastic_strain_master;
  const MaterialProperty<RankTwoTensor> &_inelastic_strain_slave;
  const MaterialProperty<RankTwoTensor> &_inelastic_strain_master_old;
  const MaterialProperty<RankTwoTensor> &_inelastic_strain_slave_old;
  ///@}
  MaterialProperty<Real> &_avg_mises_stress;
  const MaterialProperty<Real> &_avg_mises_stress_old;

  // const MaterialProperty<Real> &_avg_mises_stress_rate;
  MaterialProperty<Real> &_avg_hyd_stress;
  const MaterialProperty<Real> &_avg_hyd_stress_old;
  // const MaterialProperty<Real> &_avg_hyd_stress_rate;
  // const MaterialProperty<Real> &_avg_eq_strain;
  MaterialProperty<Real> &_avg_eq_strain_rate;
  const MaterialProperty<Real> &_avg_eq_strain_rate_old;
  MaterialProperty<Real> &_accumulated_eq_strain;
  const MaterialProperty<Real> &_accumulated_eq_strain_old;
  MaterialProperty<Real> &_interface_triaxiality;
  /* FAILURE RELATED PROPERTIES */
  MaterialProperty<bool> &_elem_failed;
  const MaterialProperty<bool> &_elem_failed_old;
  const Real _D_thr;
  const Real _max_allowed_opening_traction;
  const Real _max_allowed_damage_rate; /*maximum allowed damage rate*/
  const Real _min_allowed_residual_life;
  const Real _traction_mean_decay_time_factor;
  const Real _min_allowed_residual_stiffness;
  MaterialProperty<Real> &_time_at_failure;
  const MaterialProperty<Real> &_time_at_failure_old;
  MaterialProperty<RealVectorValue> &_traction_at_failure;
  const MaterialProperty<RealVectorValue> &_traction_at_failure_old;
  // MaterialProperty<bool> &_decay_exhausted;
  // const MaterialProperty<bool> &_decay_exhausted_old;

  MaterialProperty<RealVectorValue> &_du_at_failure;
  const MaterialProperty<RealVectorValue> &_du_at_failure_old;
  MaterialProperty<RealVectorValue> &_K_at_failure;
  const MaterialProperty<RealVectorValue> &_K_at_failure_old;

  // PETSC STUFF
  SNES _q_snes;        /* nonlinear solver context */
  KSP _q_ksp;          /* linear solver context */
  PC _q_pc;            /* preconditioner context */
  Vec _q_x, _q_r;      /* solution, residual vectors */
  Mat _q_J;            /* Jacobian matrix */
  ApplicationCtx _ctx; /* user-defined context */
  SNESConvergedReason _q_reason;
  SNESLineSearch _q_linesearch;
  SNESType _q_snestype;
  PetscErrorCode _q_ierr;
  PetscScalar *xx;
  MatScalar *mm;

  static PetscErrorCode FormJacobian1(SNES, Vec, Mat, Mat, void *);
  static PetscErrorCode FormFunction1(SNES, Vec, Vec, void *);

  void writeDataToFile() const;
  void setNewton(const bool &linesearch_on);

  const Real _mysnes_abs_tol;
  const Real _mysnes_rel_tol;
  const Real _mysnes_step_tol;
  const Real _mysnes_max_iteration;
  const bool _use_substep;
  const unsigned int _max_substep_cuts;
  const bool _force_substep;
  const bool _use_LM;
  const int _n_equation;
  const bool _give_up_qp;
  bool substepFun(nlFunBase::io_maps_type &x_sol_real, bool &fail_while_substep,
                  GBCavitationNLSystem &GBNLsystem);
  void prepareSolverContext(const Real &substep_dt,
                            nlFunBase::io_maps_type &x_old,
                            nlFunBase::io_maps_type &NL_params,
                            GBCavitationNLSystem &GBNLsystem);
  nlFunBase::io_maps_type prepareParamsSubstep(const Real &time_last_solution,
                                               const Real &current_dt);
  void initSnesGuess(const nlFunBase::io_maps_type &x_Real,
                     const nlFunBase::io_maps_type &x_Real_old,
                     const GBCavitationNLSystem &GBNLsystem);
  void updateLastSolution(nlFunBase::io_maps_type &x_last_solution) const;
  bool checkCavitationConvergence(const nlFunBase::io_maps_type &NL_params,
                                  const nlFunBase::io_maps_type &x_old,
                                  const Real &dt_local,
                                  nlFunBase::io_maps_type &x_sol_real,
                                  const GBCavitationNLSystem &GBNLsystem);
  nlFunBase::io_maps_type copyNLsolutionToMap() const;

  void updatedStateVarFromRealSolution(const nlFunBase::io_maps_type &x_Real,
                                       const Real &t_curr);

  nlFunBase::io_maps_type
  getRealSolutionFromNLSolution(const nlFunBase::io_maps_type &NL_params,
                                const nlFunBase::io_maps_type &x_old,
                                const Real &dt_local,
                                const GBCavitationNLSystem &GBNLsystem) const;

  nlFunBase::io_maps_type xoldFromOld() const;

  void decoupeldShearTraction(const Real &dt);
  void update_Dtn_dUN(const nlFunBase::io_maps_type &NL_params,
                      const nlFunBase::io_maps_type &x_old,
                      const Real &dt_local,
                      const GBCavitationNLSystem &GBNLsystem) {
    nlFunBase::io_maps_type xNL = copyNLsolutionToMap();

    DenseMatrix<Real> Jac =
        GBNLsystem.computeSystemVarJacobian(xNL, NL_params, x_old, dt_local,
                                            /*wrt_xNL=*/false);

    const std::map<std::string, DenseVector<Real>> dx_dparam =
        GBNLsystem.computeSystemParamGradient(xNL, NL_params, x_old, dt_local);
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++) {
        Jac(i, j) = -Jac(i, j);
        if (i == j)
          Jac(i, j) += 1;
      }
    DenseMatrix<Real> J_temp(3, 3);
    DenseVector<Real> dx_dy(3, 0);

    J_temp = Jac;
    DenseVector<Real> deq_dui = dx_dparam.find("un_dot")->second;
    deq_dui *= (1. / dt_local);
    J_temp.lu_solve(deq_dui, dx_dy);
    // for (unsigned int i = 0; i < 3; i++)
    _dtraction_djump[_qp](0, 0) = dx_dy(2);
    _dtraction_djump[_qp](0, 1) = 0;
    _dtraction_djump[_qp](0, 2) = 0;
  }

  void tractionDeacy(const Real &time_at_failure, const Real &a_at_failure,
                     const Real &b_at_failure, const Real &D_at_failure,
                     const Real &D_dot_at_failure, const Real &residual_life,
                     const RealVectorValue &traction_at_failure,
                     const RealVectorValue &du_at_failure,
                     const RealVectorValue &K_at_failure) {
    _elem_failed[_qp] = true;
    _time_at_failure[_qp] = time_at_failure;
    _a[_qp] = a_at_failure;
    _b[_qp] = b_at_failure;
    _D[_qp] = D_at_failure;
    _D_dot[_qp] = D_dot_at_failure;
    _residual_life[_qp] = residual_life;

    for (unsigned int i = 0; i < 3; i++) {
      _traction_at_failure[_qp](i) = traction_at_failure(i);
      _du_at_failure[_qp](i) = du_at_failure(i);
      _K_at_failure[_qp](i) = std::abs(K_at_failure(i));
    }

    Real dt_from_failure = _t - _time_at_failure[_qp];

    /*COMPUTING DECAY FACTOR*/
    Real decay_factor =
        std::exp(-(dt_from_failure) / (0.5 * _residual_life[_qp] *
                                       _traction_mean_decay_time_factor));

    std::vector<Real> K(3, 0);
    for (unsigned int i = 0; i < 3; i++)
      K[i] = std::max(_K_at_failure[_qp](i) * decay_factor,
                      _min_allowed_residual_stiffness);

    /*COMPUTING TRACTION*/
    for (unsigned int i = 0; i < 3; i++) {

      _traction[_qp](i) =
          (_displacement_jump[_qp](i) - _du_at_failure[_qp](i)) * K[i] +
          _traction_at_failure[_qp](i) * decay_factor;
    }

    /*COMPUTING DERIVATIVES*/
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        _dtraction_djump[_qp](i, j) = 0;
        if (i == j)
          _dtraction_djump[_qp](i, j) = K[i];
      }
    }
  }

  void copyOldSolution() {
    _a[_qp] = _a_old[_qp];
    _b[_qp] = _b_old[_qp];
    _residual_life[_qp] = _residual_life_old[_qp];

    for (unsigned int i = 0; i < 3; i++) {
      _traction[_qp](i) = _traction_old[_qp](i);
      _K_at_failure[_qp] = _K_at_failure_old[_qp];
    }
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        _dtraction_djump[_qp](i, j) = 0;

    _D[_qp] = _D_old[_qp];
    _elem_failed[_qp] = false;
    // _decay_exhausted[_qp] = false;
    if (_D[_qp] >= _D_thr) {
      Moose::out << "element " << _current_elem->id() << " qp " << _qp
                 << " exceeded maximum Damage, marked "
                    "as failed\n";
      _elem_failed[_qp] = true;
    }

    if (_traction[_qp](0) >= _max_allowed_opening_traction) {
      _elem_failed[_qp] = true;
      Moose::out << "element " << _current_elem->id() << " qp " << _qp
                 << " exceeded maximum Traction, "
                    "marked as failed\n";
    }

    for (unsigned int i = 0; i < 3; i++) {
      _traction_at_failure[_qp](i) = _traction[_qp](i);
      _du_at_failure[_qp](i) = _displacement_jump[_qp](i);
    }
    _time_at_failure[_qp] = _time_at_failure_old[_qp];

    _nucleation_above_threshold[_qp] = _nucleation_above_threshold_old[_qp];
    _VLdot[_qp] = _VLdot_old[_qp];
    _VHdot[_qp] = _VHdot_old[_qp];
  }
};
