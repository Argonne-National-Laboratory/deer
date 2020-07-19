//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBCavitation.h"
#include "NLSystem.h"
registerMooseObject("DeerApp", GBCavitation);

InputParameters GBCavitation::validParams() {
  InputParameters params = CZMMaterialBasePK1::validParams();
  params.addClassDescription("Sham Needleman grain boundary cavitation model. "
                             "Default parameters are for Grade91 at 600C.");
  params.addParam<Real>("a0", 5e-5, "Initial cavity radius");
  params.addParam<Real>("b0", 6e-2, "Initial cavity half spacing");
  params.addParam<Real>("FN_NI", 2e4, "Normalized nucleation rate constant");
  params.addParam<Real>("Nmax_NI", 1e3, "Normalized maximum cavity density");
  params.addParam<Real>(
      "sigma_0", 200.,
      "Traction normalization parameter for cavity nculeation");
  params.addParam<Real>("psi_degree", 75., "Cavity half-tip angle");
  params.addParam<Real>("beta", 2., "Cavity nucleation exponent");
  params.addParam<Real>("E_GB", 150e3, "Grain boundary opening stiffness");
  params.addParam<Real>(
      "E_penalty_minus_thickenss_over_2", 5,
      "Element co-penetration penatly at jump = -thickness/2");
  params.addParam<Real>("E_penalty_minus_thickenss", 10,
                        "Element co-penetration penatly at jump = -thickness");
  params.addParam<Real>("G_GB", 57.69e3, "Grain boundary shear stiffness");
  params.addParam<Real>("D_GB", 1e-15, "Grain boundary diffusivity");
  params.addParam<Real>("eta_sliding", 1e6,
                        "Grain boundary sliding viscosisty");
  params.addParam<Real>("interface_thickness", 0.0113842,
                        "The interface thickness");
  params.addParam<Real>("n", 5., "Creep rate exponent");
  params.addParam<Real>(
      "theta", 0.,
      "Parameter for the time integration scheme: 0 implicit (default), "
      " 0.5 cranck-nicolson, 1 explicit, other numbers are mixed scheme");
  params.addParam<bool>("nucleation_on", true, "Turns on cavity nucleation");
  params.addParam<bool>("growth_on", true, "Turns on cavity growth");
  params.addParam<bool>("use_triaxial_growth", true,
                        "if true(default) turns on triaxial cavity growth");
  params.addParam<bool>("use_old_bulk_property", true,
                        "If true (default) use old bulk material property.");
  params.addParam<bool>("scale_variables", true,
                        "If true (default) scale cavitation equations to "
                        "better condition the NL system.");
  params.addParam<unsigned int>("max_time_cut", 5,
                                "the maximum number of times _dt is divided "
                                "by 2 (default 5, e.g. dt_min = _dt/(2^5))");
  params.addParam<unsigned int>("max_nonlinear_iter", 10,
                                "the maximum number of nonlinear iterations "
                                "(default 10) for each substep.");
  params.addParam<Real>(
      "nl_residual_abs_tol", 1e-10,
      "The solver uses the max(abs(Residual)) as stopping criteria. 1e-10 is "
      "the default value.");
  params.addParam<Real>("D_failure", 0.95, "Damage at failure");
  params.addParam<Real>(
      "minimum_allowed_residual_life", 1e2,
      "When the interface residual life is below this value the element is "
      "marked as failed. Residual life is comptued as "
      "(1-current_damage)/damage_rate.");
  params.addParam<Real>("maximum_allowed_opening_traction", 1.4e4,
                        "When the opening traction exceeds this value the "
                        "lement is marked as failed.");
  params.addParam<Real>("minimum_allowed_stiffness", 1,
                        "The minimum allowed stiffness to prevent floating "
                        "domains during traction decay");
  params.addParam<bool>("force_substep", false, "force substep");
  return params;
}

GBCavitation::GBCavitation(const InputParameters &parameters)
    : CZMMaterialBasePK1(parameters),
      _use_old_bulk_property(getParam<bool>("use_old_bulk_property")),
      _stress_master(getMaterialPropertyByName<RankTwoTensor>("stress")),
      _stress_slave(getNeighborMaterialPropertyByName<RankTwoTensor>("stress")),
      _inelastic_strain_master(
          getMaterialPropertyByName<RankTwoTensor>("inelastic_strain")),
      _inelastic_strain_slave(
          getNeighborMaterialPropertyByName<RankTwoTensor>("inelastic_strain")),
      _inelastic_strain_master_old(
          getMaterialPropertyOld<RankTwoTensor>("inelastic_strain")),
      _inelastic_strain_slave_old(
          getNeighborMaterialPropertyOld<RankTwoTensor>("inelastic_strain")),
      _stress_vm(declareProperty<Real>("stress_vm_interface")),
      _stress_vm_old(getMaterialPropertyOldByName<Real>("stress_vm_interface")),
      _stress_H(declareProperty<Real>("stress_H_interface")),
      _stress_H_old(getMaterialPropertyOldByName<Real>("stress_H_interface")),
      _strain_rate_eq(declareProperty<Real>("strain_rate_eq_interface")),
      _strain_rate_eq_old(
          getMaterialPropertyOldByName<Real>("strain_rate_eq_interface")),
      _strain_eq(declareProperty<Real>("strain_eq_interface")),
      _strain_eq_old(getMaterialPropertyOldByName<Real>("strain_eq_interface")),
      _a(declareProperty<Real>("average_cavity_radii")),
      _a_old(getMaterialPropertyOldByName<Real>("average_cavity_radii")),
      _b(declareProperty<Real>("average_cavity_spacing")),
      _b_old(getMaterialPropertyOldByName<Real>("average_cavity_spacing")),
      _nucleation_is_active(declareProperty<int>("nucleation_is_active")),
      _nucleation_is_active_old(
          getMaterialPropertyOldByName<int>("nucleation_is_active")),
      _D(declareProperty<Real>("interface_damage")),
      _D_rate(declareProperty<Real>("interface_damage_rate")),
      _element_failed(declareProperty<int>("element_failed")),
      _element_failed_old(getMaterialPropertyOldByName<int>("element_failed")),
      _time_at_failure(declareProperty<Real>("time_at_failure")),
      _time_at_failure_old(
          getMaterialPropertyOldByName<Real>("time_at_failure")),
      _traction_at_failure(
          declareProperty<RealVectorValue>("traction_at_failure")),
      _traction_at_failure_old(
          getMaterialPropertyOldByName<RealVectorValue>("traction_at_failure")),
      _jump_at_failure(declareProperty<RealVectorValue>("jump_at_failure")),
      _jump_at_failure_old(
          getMaterialPropertyOldByName<RealVectorValue>("jump_at_failure")),
      _residual_life(declareProperty<Real>("residual_life")),
      _residual_life_old(getMaterialPropertyOldByName<Real>("residual_life")),
      _a0(getParam<Real>("a0")), _b0(getParam<Real>("b0")),
      _NI(1. / (M_PI * _b0 * _b0)), _FN_NI(getParam<Real>("FN_NI")),
      _FN(_FN_NI * _NI), _Nmax_NI(getParam<Real>("Nmax_NI")),
      _b_sat(1. / std::sqrt(M_PI * _Nmax_NI * _NI)),
      _S0(getParam<Real>("sigma_0")), _beta(getParam<Real>("beta")),
      _psi_degree(getParam<Real>("psi_degree")),
      _h(ShamNeedlemann::h_psi(_psi_degree * M_PI / 180.)),
      _E_GB(getParam<Real>("E_GB")), _G_GB(getParam<Real>("G_GB")),
      _D_GB(getParam<Real>("D_GB")),
      _eta_sliding(getParam<Real>("eta_sliding")),
      _thickness(getParam<Real>("interface_thickness")),
      _n(getParam<Real>("n")), _theta(getParam<Real>("theta")),
      _E_penalty_minus_thickenss_over_2(
          getParam<Real>("E_penalty_minus_thickenss_over_2")),
      _E_penalty_minus_thickenss(getParam<Real>("E_penalty_minus_thickenss")),
      _nucleation_on(getParam<bool>("nucleation_on")),
      _growth_on(getParam<bool>("growth_on")),
      _use_triaxial_growth(getParam<bool>("use_triaxial_growth")),
      _D_failure(getParam<Real>("D_failure")),
      _minimum_allowed_residual_life(
          getParam<Real>("minimum_allowed_residual_life")),
      _maximum_allowed_opening_traction(
          getParam<Real>("maximum_allowed_opening_traction")),
      _minimum_allowed_stiffness(getParam<Real>("minimum_allowed_stiffness")),
      _scale_variables(getParam<bool>("scale_variables")),
      _max_time_cut(getParam<unsigned int>("max_time_cut") + 1),
      _max_nonlinear_iter(getParam<unsigned int>("max_nonlinear_iter") + 1),
      _nl_residual_abs_tol(getParam<Real>("nl_residual_abs_tol")),
      _force_substep(getParam<bool>("force_substep"))

{}

void GBCavitation::computeTractionIncrementAndDerivatives() {

  if (_element_failed_old[_qp] != 0) {
    _element_failed[_qp] = 1;
    _time_at_failure[_qp] = _time_at_failure_old[_qp];
    _traction_at_failure[_qp] = _traction_at_failure_old[_qp];
    _jump_at_failure[_qp] = _jump_at_failure_old[_qp];
    _residual_life[_qp] = _residual_life_old[_qp];
    _a[_qp] = _a_old[_qp];
    _b[_qp] = _b_old[_qp];
    _D[_qp] = _a[_qp] / _b[_qp];

    tractionDecay();
  } else {
    computeAverageBulkPorperties();

    uint neq = 5;
    uint nlm = 3;
    const uint nsys = neq + nlm;

    /// set system paramters
    std::vector<string> pname, rate_pname;
    vecD pvalue, rate_pvalue;
    initNLSystemParamter(pname, pvalue, rate_pname, rate_pvalue);
    NLSystemParameters sysparams(pname, pvalue, rate_pname, rate_pvalue);

    NLVar a_var(0, "a", _a_old[_qp], _a_old[_qp]);
    NLVar b_var(1, "b", _b_old[_qp], _b_old[_qp]);
    NLVar TN_var(2, "Tn", _traction_old[_qp](0), _traction_old[_qp](0));
    NLVar TS1_var(3, "Ts1", _traction_old[_qp](1), _traction_old[_qp](1));
    NLVar TS2_var(4, "Ts2", _traction_old[_qp](2), _traction_old[_qp](2));

    /// set system up variables
    if (_scale_variables) {
      a_var.setScaleFactor(1e-3);
      b_var.setScaleFactor(1e-3);
      TN_var.setScaleFactor(1e3);
      TS1_var.setScaleFactor(1);
      TS2_var.setScaleFactor(1);
    }

    std::vector<NLVar *> myvars{&a_var, &b_var, &TN_var, &TS1_var, &TS2_var};
    NLSystemVars sysvars(myvars);

    /// set up precalculator
    ShamNeedlemann::V_dot vdotfun(&sysvars, &sysparams, {"vdot"},
                                  _use_triaxial_growth);

    /// set up equations
    ShamNeedlemann::a_res a_eq(0, sysvars, sysparams, vdotfun, _theta,
                               _growth_on);
    ShamNeedlemann::b_res b_eq(1, sysvars, sysparams, vdotfun, _theta,
                               _nucleation_on);
    ShamNeedlemann::TN_res Tn_eq(2, sysvars, sysparams, vdotfun, _theta);
    ShamNeedlemann::TS_res Ts1_eq(3, sysvars, sysparams, vdotfun, _theta, 1);
    ShamNeedlemann::TS_res Ts2_eq(4, sysvars, sysparams, vdotfun, _theta, 2);
    std::vector<Equation *> my_eqs = {&a_eq, &b_eq, &Tn_eq, &Ts1_eq, &Ts2_eq};

    /// set up lagrange multiplier equations
    ShamNeedlemann::a_lt_b c0(0, sysvars, sysparams, nsys);
    ShamNeedlemann::a_gt_a0 c1(1, sysvars, sysparams, nsys);
    ShamNeedlemann::b_lt_b_old c2(2, sysvars, sysparams, nsys);
    std::vector<const InequalityConstraint *> my_lms = {&c0, &c1, &c2};
    vecD l0(nlm);

    /// set up the constrained nonlinear system
    NLSystem mysys(&sysvars, &sysparams, my_eqs, my_lms, &vdotfun);

    vecD R = mysys.assembleR(l0);
    matrixD J = mysys.assembleJ(l0);

    /// setup the newton solver
    ShamNeedlemann::Solver newtonSolver(&mysys, &sysvars, _nl_residual_abs_tol,
                                        _max_nonlinear_iter,
                                        miconossmath::normtype::INF);

    bool converged = false;
    bool custom_interruption = false;
    Real dt_effective = 0;

    // solve the system
    matrixD deq_dparam;
    int ierr = newtonSolver.solveSubstep(
        l0, J, converged, &sysparams, {"udot_N", "udot_S1", "udot_S2"},
        deq_dparam, custom_interruption, dt_effective, _max_time_cut,
        _scale_variables, _force_substep);

    if (ierr != 0 || !converged) {
      // if we didn't converge request a global cutback
      throw MooseException(
          "GB Cavitation traction update failed, requesting global cutback");
    }
    // if we converged to a solution
    else {
      // update state var values
      _a[_qp] = a_var.getValue();
      _b[_qp] = b_var.getValue();
      _D[_qp] = _a[_qp] / _b[_qp];
      _element_failed[_qp] = sysparams.getValue("element_failed");
      _residual_life[_qp] = sysparams.getValue("residual_life");
      _time_at_failure[_qp] = _t - _dt + dt_effective;

      if (dt_effective < _dt) {

        _traction_at_failure[_qp](0) = TN_var.getValue();
        _traction_at_failure[_qp](1) = TS1_var.getValue();
        _traction_at_failure[_qp](2) = TS2_var.getValue();

        for (uint i = 0; i < 3; i++)
          if (!std::isfinite(_traction_at_failure[_qp](i)))
            mooseError("GBCavitation failed while substepping but "
                       "_traction_at_failure is "
                       "not finite");

        if (!std::isfinite(_a[_qp]))
          mooseError("GBCavitation failed while substepping but  _a[_qp] is "
                     "not finite: " +
                     std::to_string(_a[_qp]));
        if (!std::isfinite(_b[_qp]))
          mooseError("GBCavitation failed while substepping but _b[_qp]  is "
                     "not finite: " +
                     std::to_string(_b[_qp]));

        _jump_at_failure[_qp] =
            _displacement_jump_old[_qp] +
            _displacement_jump_inc[_qp] / _dt * dt_effective;

        tractionDecay();

      } else {
        _traction[_qp](0) = TN_var.getValue();
        _traction[_qp](1) = TS1_var.getValue();
        _traction[_qp](2) = TS2_var.getValue();

        for (uint i = 0; i < 3; i++)
          if (!std::isfinite(_traction[_qp](i)))
            mooseError("GBCavitation converged but traction is not finite");

        if (!std::isfinite(_a[_qp]))
          mooseError("GBCavitation converged but _a[_qp] is not finite: " +
                     std::to_string(_a[_qp]));
        if (!std::isfinite(_b[_qp]))
          mooseError("GBCavitation converged but _b[_qp] is not finite: " +
                     std::to_string(_b[_qp]));

        // update related history variables
        if (_nucleation_is_active_old[_qp])
          _nucleation_is_active[_qp] = _nucleation_is_active_old[_qp];
        else
          _nucleation_is_active[_qp] = b_eq.nucleationAboveThreshold(true);

        // record failure state

        _traction_at_failure[_qp] = _traction[_qp];
        _jump_at_failure[_qp] = _displacement_jump[_qp];

        // compute the traction increment
        _traction_inc[_qp] = _traction[_qp] - _traction_old[_qp];
        // copy derivatives in the proper container
        for (uint i = 0; i < 3; i++)
          for (uint j = 0; j < 3; j++)
            _dtraction_djump[_qp](i, j) = deq_dparam[i][j + 2];
      }
    }

    /// check for nans in traction and traction derivatives
    for (uint i = 0; i < 3; i++)
      if (!std::isfinite(_traction[_qp](i)))
        mooseError("GBCavitation _traction[_qp] " + std::to_string(i) +
                   " is not finite: " + std::to_string(_traction[_qp](i)));

    if (!std::isfinite(_a[_qp]))
      mooseError("GBCavitation _a[_qp] is not finite: " +
                 std::to_string(_a[_qp]));
    if (!std::isfinite(_b[_qp]))
      mooseError("GBCavitation _b[_qp] is not finite: " +
                 std::to_string(_b[_qp]));

    for (uint i = 0; i < 3; i++)
      for (uint j = 0; j < 3; j++)
        if (!std::isfinite(_dtraction_djump[_qp](i, j)))
          mooseError("GBCavitation _dtraction_djump[_qp](" + std::to_string(i) +
                     "," + std::to_string(j) + ") is not finite: " +
                     std::to_string(_dtraction_djump[_qp](i, j)));
  }
}

void GBCavitation::tractionDecay() {
  const Real decay_factor =
      std::exp((_time_at_failure[_qp] - _t) / (_residual_life[_qp] / 5.));
  for (int i = 0; i < 3; i++) {
    Real C = std::max(
        std::abs(_traction_at_failure[_qp](i) / _jump_at_failure[_qp](i)) *
            decay_factor,
        _minimum_allowed_stiffness);
    _traction[_qp](i) =
        (_displacement_jump[_qp](i) - _jump_at_failure[_qp](i)) * C +
        _traction_at_failure[_qp](i) * decay_factor;

    for (int j = 0; j < 3; j++) {
      if (i == j)
        _dtraction_djump[_qp](i, j) = C;
      else
        _dtraction_djump[_qp](i, j) = 0;
    }
  }
  _traction_inc[_qp] = _traction[_qp] - _traction_old[_qp];
}

void GBCavitation::initQpStatefulProperties() {
  CZMMaterialBasePK1::initQpStatefulProperties();
  _stress_vm[_qp] = 0;
  _stress_H[_qp] = 0;
  _strain_rate_eq[_qp] = 0;
  _a[_qp] = _a0;
  _b[_qp] = _b0;
  _nucleation_is_active[_qp] = 0;
  _element_failed[_qp] = 0;
  _time_at_failure[_qp] = 0;
  _traction_at_failure[_qp] = 0;
  _jump_at_failure[_qp] = 0;
  _residual_life[_qp] = 1e6;
}

void GBCavitation::computeAverageBulkPorperties() {
  // compute average Von Mises Stress;
  RankTwoTensor dev_stress = _stress_master[_qp].deviatoric();
  _stress_vm[_qp] = std::sqrt(1.5 * dev_stress.doubleContraction(dev_stress));
  dev_stress = _stress_slave[_qp].deviatoric();
  _stress_vm[_qp] += std::sqrt(1.5 * dev_stress.doubleContraction(dev_stress));
  _stress_vm[_qp] /= 2.;

  // compute average hydrostatic stress;
  _stress_H[_qp] =
      (_stress_master[_qp].trace() + _stress_slave[_qp].trace()) / 6.;

  // compute equivalent inelastic strain rate;
  RankTwoTensor strain_rate =
      (_inelastic_strain_master[_qp] - _inelastic_strain_master_old[_qp]) / _dt;
  _strain_rate_eq[_qp] =
      std::sqrt(2.0 / 3.0 * strain_rate.doubleContraction(strain_rate));
  strain_rate =
      (_inelastic_strain_slave[_qp] - _inelastic_strain_slave_old[_qp]) / _dt;
  _strain_rate_eq[_qp] +=
      std::sqrt(2.0 / 3.0 * strain_rate.doubleContraction(strain_rate));
  _strain_rate_eq[_qp] /= 2.;

  _strain_eq[_qp] = _strain_eq_old[_qp] + _strain_rate_eq[_qp] * _dt;
}

void GBCavitation::initNLSystemParamter(std::vector<std::string> &pname,
                                        vecD &pvalue,
                                        std::vector<std::string> &rate_pname,
                                        vecD &rate_pvalue) {
  pname = {"dt",
           "dt_accum",
           "max_ab",
           "FN",
           "FN_NI",
           "a0",
           "b_sat",
           "S0",
           "beta",
           "h",
           "E",
           "G",
           "D",
           "nucleation_is_active",
           "n",
           "eta_sliding",
           "thickness",
           "Svm",
           "Sh",
           "e",
           "max_damage",
           "minimum_allowed_residual_life",
           "maximum_allowed_opening_traction",
           "residual_life",
           "element_failed",
           "uN_old"};
  pvalue = {_dt,
            _dt,
            1. - 1e-3,
            _FN,
            _FN_NI,
            _a0,
            _b_sat,
            _S0,
            _beta,
            ShamNeedlemann::h_psi(_psi_degree * M_PI / 180.),
            _E_GB,
            _G_GB,
            _D_GB,
            (double)_nucleation_is_active[_qp],
            _n,
            _eta_sliding,
            _thickness,
            _use_old_bulk_property ? _stress_vm_old[_qp] : _stress_vm[_qp],
            _use_old_bulk_property ? _stress_H_old[_qp] : _stress_H[_qp],
            _use_old_bulk_property ? _strain_eq_old[_qp] : _strain_eq[_qp],
            _D_failure,
            _minimum_allowed_residual_life,
            _maximum_allowed_opening_traction,
            1e6,
            0.,
            _displacement_jump_old[_qp](0)};

  rate_pname = {"udot_N", "udot_S1", "udot_S2", "edot"};
  rate_pvalue = {_displacement_jump_inc[_qp](0) / _dt,
                 _displacement_jump_inc[_qp](1) / _dt,
                 _displacement_jump_inc[_qp](2) / _dt,
                 _use_old_bulk_property ? _strain_rate_eq_old[_qp]
                                        : _strain_rate_eq[_qp]};
}

// void GBCavitation::decoupeldShearTraction(const Real &dt) {
//   Real a = _a_old[_qp];
//   Real b = _b_old[_qp];
//   Real S = _thickness / (_G_GB * (1. - (a / b)));
//   Real eta = _eta_sliding;
//   if (a / b > 0.5)
//     eta = 2 * _eta_sliding * (1. - a / b);
//   _traction[_qp](1) =
//       eta * _displacement_jump_inc[_qp](1) / _dt +
//       std::exp(-dt / (eta * S)) *
//           (_traction_old[_qp](1) - eta * _displacement_jump_inc[_qp](1) /
//           _dt);
//
//   _dtraction_djump[_qp](1, 0) = 0;
//   _dtraction_djump[_qp](1, 1) = eta / dt - std::exp(-dt / (eta * S)) * eta
//   / dt; _dtraction_djump[_qp](1, 2) = 0;
//
//   _traction[_qp](2) =
//       eta * _displacement_jump_inc[_qp](2) / _dt +
//       std::exp(-dt / (eta * S)) *
//           (_traction_old[_qp](2) - eta * _displacement_jump_inc[_qp](2) /
//           _dt);
//
//   _dtraction_djump[_qp](2, 0) = 0;
//   _dtraction_djump[_qp](2, 1) = 0;
//   _dtraction_djump[_qp](2, 2) = eta / dt - std::exp(-dt / (eta * S)) * eta
//   / dt;
// }
