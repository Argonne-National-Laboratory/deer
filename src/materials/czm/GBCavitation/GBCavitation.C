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

InputParameters
GBCavitation::validParams()
{
  InputParameters params = CZMComputeLocalTractionIncrementalBase::validParams();
  params.addClassDescription("Sham Needleman grain boundary cavitation model. "
                             "Default parameters are for Grade91 at 600C.");
  params.addParam<Real>("a0", 5e-5, "Initial cavity radius, a0");
  params.addParam<Real>("b0", 6e-2, "Initial cavity half spacing, b0");
  params.addParam<Real>("FN_NI", 2e4, "Normalized nucleation rate constant");
  params.addParam<Real>("Nmax_NI", 1e3, "Normalized maximum cavity density constant");
  params.addParam<Real>("sigma_0", 200., "Traction normalization parameter for cavity nucleation");
  params.addParam<Real>("psi_degree", 75., "Cavity half-tip angle");
  params.addParam<Real>("beta_exponent", 2., "Cavity nucleation exponent");
  params.addParam<Real>("E_GB", 150e3, "Grain boundary opening modulus (stress)");
  params.addParam<Real>("G_GB", 58.3657588e3, "Grain boundary sliding modulus (stress)");
  params.addParam<Real>(
      "E_penalty_minus_thickenss", 10, "Element co-penetration penatly at jump = -thickness");
  params.addParam<Real>("E_penalty_after_failure_minus_thickenss",
                        1e6,
                        "Element co-penetration penatly at jump = -thickness");
  params.addParam<Real>("D_GB", 1e-15, "Grain boundary diffusivity");
  params.addParam<Real>("eta_sliding", 1e6, "Grain boundary sliding viscosisty");
  params.addParam<Real>("interface_thickness", 0.0113842, "The interface thickness");
  params.addParam<Real>("interface_thickness_after_failure",
                        0.00113842,
                        "The interface thickness for the failure model");
  params.addParam<Real>("n_exponent", 5., "Creep rate exponent");
  params.addParam<Real>("theta",
                        0.,
                        "Parameter for the time integration scheme: 0 implicit (default), "
                        " 0.5 cranck-nicolson, 1 explicit, other numbers are mixed scheme");
  params.addParam<bool>("nucleation_on", true, "Turns on cavity nucleation");
  params.addParam<bool>("growth_on", true, "Turns on cavity growth");
  params.addParam<bool>(
      "use_triaxial_growth", true, "if true(default) turns on triaxial cavity growth");
  params.addParam<unsigned int>("vdot_method", 0, "0 use vdotL, 1 use vdotH, 2 use both");
  params.addParam<bool>("use_old_bulk_property", false, "If true use old bulk material property.");
  params.addParam<unsigned int>("max_time_cut",
                                5,
                                "the maximum number of times _dt is divided "
                                "by 2 (default 5, e.g. dt_min = _dt/(2^5))");
  params.addParam<unsigned int>("max_nonlinear_iter",
                                10,
                                "the maximum number of nonlinear iterations "
                                "(default 10) for each substep.");
  params.addParam<Real>("nl_residual_abs_tol",
                        1e-12,
                        "The solver uses the max(abs(Residual)) as stopping criteria. 1e-12 is "
                        "the default value.");
  params.addParam<Real>("D_failure", 0.95, "Damage at failure, a/b");
  params.addParam<Real>("minimum_allowed_residual_life",
                        1e2,
                        "When the interface residual life is below this value the element is "
                        "marked as failed. Residual life is comptued as "
                        "(1-current_damage)/damage_rate.");
  params.addParam<Real>("maximum_allowed_opening_traction",
                        1.4e4,
                        "When the opening traction exceeds this value the "
                        "lement is marked as failed.");
  params.addParam<Real>("minimum_allowed_stiffness",
                        1,
                        "The minimum allowed stiffness to prevent floating "
                        "domains during traction decay");
  params.addParam<bool>("force_substep", false, "force substep");
  params.addParam<UserObjectName>(
      "GBCavitationBoundaryPropertyUO",
      "the user object containing material GB dependent material properties");
  params.addParam<Real>(
      "JxW_ref", 1.9e-5, "The reference qp area used to normalize the interface stiffness");
  params.addParam<bool>("scale_stiffness",
                        false,
                        "If true (default) scale the interface stiffness according to hmin");
  return params;
}

GBCavitation::GBCavitation(const InputParameters & parameters)
  : CZMComputeLocalTractionIncrementalBase(parameters),
    _ndisp(_mesh.dimension()),
    _use_old_bulk_property(getParam<bool>("use_old_bulk_property")),
    _stress_master(getMaterialPropertyByName<RankTwoTensor>("cauchy_stress")),
    _stress_slave(getNeighborMaterialPropertyByName<RankTwoTensor>("cauchy_stress")),
    _inelastic_strain_master(getMaterialPropertyByName<RankTwoTensor>("inelastic_strain")),
    _inelastic_strain_slave(getNeighborMaterialPropertyByName<RankTwoTensor>("inelastic_strain")),
    _inelastic_strain_master_old(getMaterialPropertyOld<RankTwoTensor>("inelastic_strain")),
    _inelastic_strain_slave_old(getNeighborMaterialPropertyOld<RankTwoTensor>("inelastic_strain")),
    _stress_vm(declareProperty<Real>("stress_vm_interface")),
    _stress_vm_old(getMaterialPropertyOldByName<Real>("stress_vm_interface")),
    _stress_H(declareProperty<Real>("stress_H_interface")),
    _stress_H_old(getMaterialPropertyOldByName<Real>("stress_H_interface")),
    _strain_rate_eq(declareProperty<Real>("strain_rate_eq_interface")),
    _strain_rate_eq_old(getMaterialPropertyOldByName<Real>("strain_rate_eq_interface")),
    _strain_eq(declareProperty<Real>("strain_eq_interface")),
    _strain_eq_old(getMaterialPropertyOldByName<Real>("strain_eq_interface")),
    _a(declareProperty<Real>("a")),
    _a_old(getMaterialPropertyOldByName<Real>("a")),
    _b(declareProperty<Real>("b")),
    _b_old(getMaterialPropertyOldByName<Real>("b")),
    _nucleation_is_active(declareProperty<int>("nucleation_is_active")),
    _nucleation_is_active_old(getMaterialPropertyOldByName<int>("nucleation_is_active")),
    _D(declareProperty<Real>("interface_damage")),
    _D_rate(declareProperty<Real>("interface_damage_rate")),
    _VLdot(declareProperty<Real>("VLdot")),
    _VL1dot(declareProperty<Real>("VL1dot")),
    _VL2dot(declareProperty<Real>("VL2dot")),
    _L(declareProperty<Real>("L")),
    _element_failed(declareProperty<int>("element_failed")),
    _element_failed_old(getMaterialPropertyOldByName<int>("element_failed")),
    _time_at_failure(declareProperty<Real>("time_at_failure")),
    _time_at_failure_old(getMaterialPropertyOldByName<Real>("time_at_failure")),
    _traction_at_failure(declareProperty<RealVectorValue>("traction_at_failure")),
    _traction_at_failure_old(getMaterialPropertyOldByName<RealVectorValue>("traction_at_failure")),
    _jump_at_failure(declareProperty<RealVectorValue>("jump_at_failure")),
    _jump_at_failure_old(getMaterialPropertyOldByName<RealVectorValue>("jump_at_failure")),
    _residual_life(declareProperty<Real>("residual_life")),
    _residual_life_old(getMaterialPropertyOldByName<Real>("residual_life")),
    _n_exponent(declareProperty<Real>("n_exponent")),
    _n_exponent_old(getMaterialPropertyOldByName<Real>("n_exponent")),
    _beta_exponent(declareProperty<Real>("beta_exponent")),
    _beta_exponent_old(getMaterialPropertyOldByName<Real>("beta_exponent")),
    _a0(declareProperty<Real>("a0")),
    _a0_old(getMaterialPropertyOldByName<Real>("a0")),
    _NI(declareProperty<Real>("NI")),
    _NI_old(getMaterialPropertyOldByName<Real>("NI")),
    _FN(declareProperty<Real>("FN")),
    _FN_old(getMaterialPropertyOldByName<Real>("FN")),
    _D_GB(declareProperty<Real>("GB_diffusivity")),
    _D_GB_old(getMaterialPropertyOldByName<Real>("GB_diffusivity")),
    _eta_sliding(declareProperty<Real>("GB_sliding_viscosity")),
    _eta_sliding_old(getMaterialPropertyOldByName<Real>("GB_sliding_viscosity")),
    _h(declareProperty<Real>("h")),
    _h_old(getMaterialPropertyOldByName<Real>("h")),
    _b_sat(declareProperty<Real>("saturation_cavity_spacing")),
    _b_sat_old(getMaterialPropertyOldByName<Real>("saturation_cavity_spacing")),
    _E_GB(declareProperty<Real>("GB_young_modulus")),
    _E_GB_old(getMaterialPropertyOldByName<Real>("GB_young_modulus")),
    _G_GB(declareProperty<Real>("GB_shear_modulus")),
    _G_GB_old(getMaterialPropertyOldByName<Real>("GB_shear_modulus")),
    _thickness(declareProperty<Real>("GB_thickness")),
    _thickness_old(getMaterialPropertyOldByName<Real>("GB_thickness")),
    _sigma_0(declareProperty<Real>("sigma_0")),
    _sigma_0_old(getMaterialPropertyOldByName<Real>("sigma_0")),
    _E_penalty_minus_thickenss(getParam<Real>("E_penalty_minus_thickenss")),
    _E_penalty_after_failure_minus_thickenss(
        getParam<Real>("E_penalty_after_failure_minus_thickenss")),
    _thickness_after_failure(getParam<Real>("interface_thickness_after_failure")),
    _theta(getParam<Real>("theta")),
    _nucleation_on(getParam<bool>("nucleation_on")),
    _growth_on(getParam<bool>("growth_on")),
    _use_triaxial_growth(getParam<bool>("use_triaxial_growth")),
    _vdot_method(getParam<unsigned int>("vdot_method")),
    _D_failure(getParam<Real>("D_failure")),
    _minimum_allowed_residual_life(getParam<Real>("minimum_allowed_residual_life")),
    _maximum_allowed_opening_traction(getParam<Real>("maximum_allowed_opening_traction")),
    _minimum_allowed_stiffness(getParam<Real>("minimum_allowed_stiffness")),
    _max_time_cut(getParam<unsigned int>("max_time_cut") + 1),
    _max_nonlinear_iter(getParam<unsigned int>("max_nonlinear_iter") + 1),
    _nl_residual_abs_tol(getParam<Real>("nl_residual_abs_tol")),
    _force_substep(getParam<bool>("force_substep")),
    _GBCavitationBoundaryPropertyUO(
        parameters.isParamSetByUser("GBCavitationBoundaryPropertyUO")
            ? &getUserObjectByName<GBCavitationBoundaryPropertyUO>(
                  getParam<UserObjectName>("GBCavitationBoundaryPropertyUO"))
            : nullptr),
    _JxW_ref(getParam<Real>("JxW_ref")),
    _scale_stiffness(getParam<bool>("scale_stiffness"))
{
  // sanity checks
  if (_E_penalty_after_failure_minus_thickenss <= 1)
    mooseError("E_penalty_after_failure_minus_thickenss must be greater or "
               "equalt to 1");

  if (_E_penalty_minus_thickenss <= 1)
    mooseError("E_penalty_minus_thickenss must be greater or equalt to 1");

  if (_GBCavitationBoundaryPropertyUO == nullptr)
  {
    if (getParam<Real>("psi_degree") < 0 || getParam<Real>("psi_degree") > 90.)
      mooseError("psi_degree must be between 0 and 90 degrees");
  }
  if (_JxW_ref <= 0)
    mooseError("the paremter h_ref must be positive");

  if (_ndisp == 1)
    mooseError("not implemented for 1D problems");
}

void
GBCavitation::updateGBDependentProperties()
{
  // copy qp dependent properties
  _a0[_qp] = _a0_old[_qp];
  _NI[_qp] = _NI_old[_qp];
  _FN[_qp] = _FN_old[_qp];
  _D_GB[_qp] = _D_GB_old[_qp];
  _eta_sliding[_qp] = _eta_sliding_old[_qp];
  _h[_qp] = _h_old[_qp];
  _b_sat[_qp] = _b_sat_old[_qp];
  _E_GB[_qp] = _E_GB_old[_qp];
  _G_GB[_qp] = _G_GB_old[_qp];
  _thickness[_qp] = _thickness_old[_qp];
  _sigma_0[_qp] = _sigma_0_old[_qp];
  _n_exponent[_qp] = _n_exponent_old[_qp];
  _beta_exponent[_qp] = _beta_exponent_old[_qp];
}

void
GBCavitation::updateFailedElementProperties()
{
  // copy qp dependent properties
  _element_failed[_qp] = 1;
  _time_at_failure[_qp] = _time_at_failure_old[_qp];
  _traction_at_failure[_qp] = _traction_at_failure_old[_qp];
  _jump_at_failure[_qp] = _jump_at_failure_old[_qp];
  _residual_life[_qp] = _residual_life_old[_qp];
  _a[_qp] = _a_old[_qp];
  _b[_qp] = _b_old[_qp];
  _D[_qp] = _a[_qp] / _b[_qp];
}

NLSystemParameters
GBCavitation::setupLinearSystemParameters()
{
  /// set system paramters
  std::vector<string> pname, rate_pname;
  vecD pvalue, rate_pvalue;
  initNLSystemParamter(pname, pvalue, rate_pname, rate_pvalue);
  NLSystemParameters sysparams(pname, pvalue, rate_pname, rate_pvalue);
  return sysparams;
}

NLSystemVars
GBCavitation::setupLinearSystemVariables()
{
  NLSystemVars sysvars({"a", "b", "Tn", "Ts1", "Ts2"});
  sysvars.setValueOld("a", _a_old[_qp]);
  sysvars.setValueOld("b", _b_old[_qp]);
  sysvars.setValueOld("Tn", _interface_traction_old[_qp](0));
  sysvars.setValueOld("Ts1", _interface_traction_old[_qp](1));
  sysvars.setValueOld("Ts2", _interface_traction_old[_qp](2));
  sysvars.setToOld();
  return sysvars;
}

void
GBCavitation::setupShamNeedlemanEquations(std::vector<Equation *> & sys_equations,
                                          NLSystemVars & sysvars,
                                          NLSystemParameters & sysparams,
                                          ShamNeedlemann::V_dot & vdotfun)
{
  sys_equations.push_back(new ShamNeedlemann::a_res(
      0, sysvars, sysparams, vdotfun, _h[_qp], _a0[_qp], _theta, _growth_on));
  sys_equations.push_back(new ShamNeedlemann::b_res(1,
                                                    sysvars,
                                                    sysparams,
                                                    vdotfun,
                                                    _FN[_qp],
                                                    _FN[_qp] / _NI[_qp],
                                                    _sigma_0[_qp],
                                                    _beta_exponent[_qp],
                                                    _b_sat[_qp],
                                                    _theta,
                                                    _nucleation_on));
  sys_equations.push_back(new ShamNeedlemann::TN_res(2,
                                                     sysvars,
                                                     sysparams,
                                                     vdotfun,
                                                     _thickness[_qp],
                                                     _E_GB[_qp],
                                                     _E_penalty_minus_thickenss,
                                                     _thickness[_qp],
                                                     _theta));
  sys_equations.push_back(new ShamNeedlemann::TS_res(
      3, sysvars, sysparams, vdotfun, 1, _thickness[_qp], _eta_sliding[_qp], _G_GB[_qp], _theta));
  sys_equations.push_back(new ShamNeedlemann::TS_res(
      4, sysvars, sysparams, vdotfun, 2, _thickness[_qp], _eta_sliding[_qp], _G_GB[_qp], _theta));
}

void
GBCavitation::setupShamNeedlemanConstraints(std::vector<const InequalityConstraint *> & my_lms,
                                            NLSystemVars & sysvars,
                                            NLSystemParameters & sysparams)
{
  my_lms.push_back(new ShamNeedlemann::a_lt_b(0, sysvars, sysparams, 8));
  my_lms.push_back(new ShamNeedlemann::a_gt_a0(1, sysvars, sysparams, 8, _a0[_qp]));
  my_lms.push_back(new ShamNeedlemann::b_lt_b_old(2, sysvars, sysparams, 8));
}

NLSystem
GBCavitation::setupNonLinearSystem(std::vector<Equation *> & sys_equations,
                                   NLSystemVars & sysvars,
                                   NLSystemParameters & sysparams,
                                   ShamNeedlemann::V_dot & vdotfun,
                                   std::vector<const InequalityConstraint *> & my_lms)
{

  /// set up the constrained nonlinear system
  return NLSystem(&sysvars, &sysparams, sys_equations, my_lms, &vdotfun);
}

ShamNeedlemann::Solver
GBCavitation::setupNewtonSolver(NLSystem & mysys, NLSystemVars & sysvars)
{

  return ShamNeedlemann::Solver(
      &mysys, &sysvars, _nl_residual_abs_tol, _max_nonlinear_iter, miconossmath::normtype::INF);
}

void
GBCavitation::updateVariablesAfterNonLinearSolution(NLSystemVars & sysvars,
                                                    NLSystemParameters & sysparams,
                                                    ShamNeedlemann::V_dot & vdotfun,
                                                    Real & dt_effective)
{
  _a[_qp] = sysvars.getValue("a");
  _b[_qp] = sysvars.getValue("b");
  _D[_qp] = _a[_qp] / _b[_qp];
  _element_failed[_qp] = sysparams.getValue("element_failed");
  _residual_life[_qp] = sysparams.getValue("residual_life");
  _time_at_failure[_qp] = _t - _dt + dt_effective;

  _VLdot[_qp] = vdotfun.getValue("vdot", true);
  _VL1dot[_qp] = vdotfun.getValue("vL1dot", true);
  _VL2dot[_qp] = vdotfun.getValue("vL2dot", true);
  _L[_qp] = vdotfun.getValue("L", true);
}

void
GBCavitation::updateIfElementFailedWhilseSubstepping(NLSystemVars & sysvars, Real & dt_effective)
{

  _traction_at_failure[_qp](0) = sysvars.getValue("Tn");
  _traction_at_failure[_qp](1) = sysvars.getValue("Ts1");
  if (_ndisp == 3)
    _traction_at_failure[_qp](2) = sysvars.getValue("Ts2");

  _jump_at_failure[_qp] = _interface_displacement_jump_old[_qp] +
                          _interface_displacement_jump_inc[_qp] / _dt * dt_effective;

  // some checks, only performed in debug or devel mode
  for (uint i = 0; i < 3; i++)
  {
    mooseAssert(std::isfinite(_traction_at_failure[_qp](i)),
                "GBCavitation failed while substepping and "
                "_traction_at_failure is "
                "not finite");
  }

  mooseAssert(std::isfinite(_a[_qp]),
              "GBCavitation failed while substepping and _a[_qp] is "
              "not finite: " +
                  std::to_string(_a[_qp]));

  mooseAssert(std::isfinite(_b[_qp]),
              "GBCavitation failed while substepping and _b[_qp]  is "
              "not finite: " +
                  std::to_string(_b[_qp]));

  tractionDecay();
}

void
GBCavitation::updateForFullStep(NLSystemVars & sysvars,
                                Real & dt_effective,
                                std::vector<Equation *> & sys_equations,
                                matrixD & deq_dparam)
{
  _interface_traction[_qp](0) = sysvars.getValue("Tn");
  _interface_traction[_qp](1) = sysvars.getValue("Ts1");
  if (_ndisp == 3)
    _interface_traction[_qp](2) = sysvars.getValue("Ts2");

  // some checks, only performed in debug or devel mode
  for (uint i = 0; i < 3; i++)
  {
    mooseAssert(std::isfinite(_interface_traction[_qp](i)),
                "GBCavitation converged but traction is not finite");
  }

  mooseAssert(std::isfinite(_a[_qp]),
              "GBCavitation converged but _a[_qp] is not finite: " + std::to_string(_a[_qp]));

  mooseAssert(std::isfinite(_b[_qp]),
              "GBCavitation converged but _b[_qp] is not finite: " + std::to_string(_b[_qp]));

  // update related history variables
  if (_nucleation_is_active_old[_qp])
    _nucleation_is_active[_qp] = _nucleation_is_active_old[_qp];
  else
    _nucleation_is_active[_qp] =
        dynamic_cast<ShamNeedlemann::b_res *>(sys_equations[1])->nucleationAboveThreshold(true);

  // record failure state
  _traction_at_failure[_qp] = _interface_traction[_qp];
  _jump_at_failure[_qp] = _interface_displacement_jump[_qp];

  // compute the traction increment
  _interface_traction_inc[_qp] = _interface_traction[_qp] - _interface_traction_old[_qp];
  // copy derivatives in the proper container
  for (uint i = 0; i < 3; i++)
    for (uint j = 0; j < 3; j++)
      _dinterface_traction_djump[_qp](i, j) = deq_dparam[i][j + 2];
}

void
GBCavitation::postSolutionDebugChecks()
{
  /// check for nans in traction and traction derivatives, only performed in
  /// debug or devel mode
  for (uint i = 0; i < 3; i++)
    mooseAssert(std::isfinite(_interface_traction[_qp](i)),
                "GBCavitation _interface_traction[_qp] " + std::to_string(i) +
                    " is not finite: " + std::to_string(_interface_traction[_qp](i)));

  mooseAssert(std::isfinite(_a[_qp]),
              "GBCavitation _a[_qp] is not finite: " + std::to_string(_a[_qp]));

  mooseAssert(std::isfinite(_b[_qp]),
              "GBCavitation _b[_qp] is not finite: " + std::to_string(_b[_qp]));

  for (uint i = 0; i < 3; i++)
  {
    for (uint j = 0; j < 3; j++)
    {
      mooseAssert(std::isfinite(_dinterface_traction_djump[_qp](i, j)),
                  "GBCavitation _dinterface_traction_djump[_qp](" + std::to_string(i) + "," +
                      std::to_string(j) +
                      ") is not finite: " + std::to_string(_dinterface_traction_djump[_qp](i, j)) +
                      " increment i = " + std::to_string(_interface_displacement_jump_inc[_qp](i)));
    }
  }
  mooseAssert(std::isfinite(_dinterface_traction_djump[_qp].det()),
              "determinant of _dinterface_traction_djump is not finite");
}

void
GBCavitation::computeInterfaceTractionIncrementAndDerivatives()
{
  if (_t == 0.)
    return;

  updateGBDependentProperties();

  if (_element_failed_old[_qp] != 0)
  {
    updateFailedElementProperties();
    tractionDecay();
  }
  else
  {
    computeAverageBulkProperties();

    NLSystemParameters sysparams = setupLinearSystemParameters();
    NLSystemVars sysvars = setupLinearSystemVariables();

    /// set up precalculator
    ShamNeedlemann::V_dot vdotfun(&sysvars,
                                  &sysparams,
                                  {"vdot", "vL1dot", "vL2dot", "vH1dot", "vH2dot", "L"},
                                  _n_exponent[_qp],
                                  _h[_qp],
                                  _D_GB[_qp],
                                  _use_triaxial_growth,
                                  _vdot_method,
                                  _nucleation_on);

    std::vector<Equation *> my_eqs;
    setupShamNeedlemanEquations(my_eqs, sysvars, sysparams, vdotfun);

    std::vector<const InequalityConstraint *> my_lms;
    setupShamNeedlemanConstraints(my_lms, sysvars, sysparams);

    NLSystem mysys = setupNonLinearSystem(my_eqs, sysvars, sysparams, vdotfun, my_lms);
    ShamNeedlemann::Solver newtonSolver = setupNewtonSolver(mysys, sysvars);

    /// declare some useful variables and actually solve the problem
    bool converged = false;
    bool custom_interruption = false;
    Real dt_effective = 0;
    matrixD J(8, vecD(8));
    vecD l0(3);

    matrixD deq_dparam;
    int ierr = newtonSolver.solveSubstep(l0,
                                         J,
                                         converged,
                                         &sysparams,
                                         {"udot_N", "udot_S1", "udot_S2"},
                                         deq_dparam,
                                         custom_interruption,
                                         dt_effective,
                                         _max_time_cut,
                                         /*auto_scale_equation = */ false,
                                         _force_substep);

    if (ierr != 0 || !converged)
    {
      // if we didn't converge request a global cutback
      throw MooseException("GB Cavitation traction update failed, requesting global cutback");
    }
    // if we converged to a solution
    else
    {
      updateVariablesAfterNonLinearSolution(sysvars, sysparams, vdotfun, dt_effective);
      if (dt_effective < _dt)
        updateIfElementFailedWhilseSubstepping(sysvars, dt_effective);
      else
        updateForFullStep(sysvars, dt_effective, my_eqs, deq_dparam);
    }
    postSolutionDebugChecks();
  }
}

void
GBCavitation::tractionDecay()
{
  Real P = 1.;
  Real dP_djump = 0.;
  if (_interface_displacement_jump[_qp](0) < 0)
  {
    const Real jump = _interface_displacement_jump[_qp](0);
    const Real & P_mt = _E_penalty_after_failure_minus_thickenss;
    const double a_parabola = (P_mt - 1.) / (_thickness_after_failure * _thickness_after_failure);
    P = a_parabola * jump * jump + 1;
    dP_djump = 2 * a_parabola * jump;
  }

  const Real decay_factor = std::exp((_time_at_failure[_qp] - _t) / (_residual_life[_qp] / 5.));
  for (unsigned int i = 0; i < _ndisp; i++)
  {
    Real C =
        std::max(std::abs(_traction_at_failure[_qp](i) / _jump_at_failure[_qp](i)) * decay_factor,
                 _minimum_allowed_stiffness);
    if (i == 0)
      C *= P;
    _interface_traction[_qp](i) =
        (_interface_displacement_jump[_qp](i) - _jump_at_failure[_qp](i) * decay_factor) * C +
        _traction_at_failure[_qp](i) * decay_factor;

    for (unsigned int j = 0; j < _ndisp; j++)
    {
      if (i == j)
      {
        _dinterface_traction_djump[_qp](i, j) = C;
        if (i == 0)
        {
          _dinterface_traction_djump[_qp](i, j) +=
              (_interface_displacement_jump[_qp](i) - _jump_at_failure[_qp](i) * decay_factor) * C /
              P * dP_djump;
        }
      }
      else
        _dinterface_traction_djump[_qp](i, j) = 0;
    }
  }

  _interface_traction_inc[_qp] = _interface_traction[_qp] - _interface_traction_old[_qp];
}

void
GBCavitation::initQpStatefulProperties()
{
  CZMComputeLocalTractionIncrementalBase::initQpStatefulProperties();
  _stress_vm[_qp] = 0;
  _stress_H[_qp] = 0;
  _strain_rate_eq[_qp] = 0;
  _nucleation_is_active[_qp] = 0;
  _element_failed[_qp] = 0;
  _time_at_failure[_qp] = 0;
  _traction_at_failure[_qp] = 0;
  _jump_at_failure[_qp] = 0;
  _residual_life[_qp] = 1e6;

  InitGBCavitationParamsAndProperties();
}

void
GBCavitation::computeAverageBulkProperties()
{
  // compute average Von Mises Stress;
  RankTwoTensor dev_stress = _stress_master[_qp].deviatoric();
  _stress_vm[_qp] = std::sqrt(1.5 * dev_stress.doubleContraction(dev_stress));
  dev_stress = _stress_slave[_qp].deviatoric();
  _stress_vm[_qp] += std::sqrt(1.5 * dev_stress.doubleContraction(dev_stress));
  _stress_vm[_qp] /= 2.;

  // compute average hydrostatic stress;
  _stress_H[_qp] = (_stress_master[_qp].trace() + _stress_slave[_qp].trace()) / 6.;

  // compute equivalent inelastic strain rate;
  RankTwoTensor strain_rate =
      (_inelastic_strain_master[_qp] - _inelastic_strain_master_old[_qp]) / _dt;
  _strain_rate_eq[_qp] = std::sqrt(2.0 / 3.0 * strain_rate.doubleContraction(strain_rate));
  strain_rate = (_inelastic_strain_slave[_qp] - _inelastic_strain_slave_old[_qp]) / _dt;
  _strain_rate_eq[_qp] += std::sqrt(2.0 / 3.0 * strain_rate.doubleContraction(strain_rate));
  _strain_rate_eq[_qp] /= 2.;

  _strain_eq[_qp] = _strain_eq_old[_qp] + _strain_rate_eq[_qp] * _dt;
}

void
GBCavitation::initNLSystemParamter(std::vector<std::string> & pname,
                                   vecD & pvalue,
                                   std::vector<std::string> & rate_pname,
                                   vecD & rate_pvalue)
{
  pname = {"dt",
           "dt_accum",
           "max_ab",
           "nucleation_is_active",
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
            (double)_nucleation_is_active[_qp],
            _use_old_bulk_property ? _stress_vm_old[_qp] : _stress_vm[_qp],
            _use_old_bulk_property ? _stress_H_old[_qp] : _stress_H[_qp],
            _use_old_bulk_property ? _strain_eq_old[_qp] : _strain_eq[_qp],
            _D_failure,
            _minimum_allowed_residual_life,
            _maximum_allowed_opening_traction,
            1e6,
            0.,
            _interface_displacement_jump_old[_qp](0)};

  rate_pname = {"udot_N", "udot_S1", "udot_S2", "edot"};
  rate_pvalue = {_interface_displacement_jump_inc[_qp](0) / _dt,
                 _interface_displacement_jump_inc[_qp](1) / _dt,
                 _ndisp == 3 ? _interface_displacement_jump_inc[_qp](2) / _dt : 0.,
                 _use_old_bulk_property ? _strain_rate_eq_old[_qp] : _strain_rate_eq[_qp]};
}

void
GBCavitation::InitGBCavitationParamsAndProperties()
{
  Real FN_NI, Nmax_NI, a0, b0, psi, D_GB, E_GB, G_GB, eta_sliding, sigma_0, thickness,
      beta_exponent, n_exponent;
  if (_GBCavitationBoundaryPropertyUO == nullptr)
    getInitPropertyValuesFromParams(FN_NI,
                                    Nmax_NI,
                                    a0,
                                    b0,
                                    psi,
                                    D_GB,
                                    E_GB,
                                    G_GB,
                                    eta_sliding,
                                    sigma_0,
                                    thickness,
                                    beta_exponent,
                                    n_exponent);
  else
    getInitPropertyValuesFromUO(FN_NI,
                                Nmax_NI,
                                a0,
                                b0,
                                psi,
                                D_GB,
                                E_GB,
                                G_GB,
                                eta_sliding,
                                sigma_0,
                                thickness,
                                beta_exponent,
                                n_exponent);

  psi *= M_PI / 180.;

  _a0[_qp] = a0;
  _a[_qp] = a0;
  _b[_qp] = b0;
  _NI[_qp] = (1. / (M_PI * b0 * b0));
  _FN[_qp] = FN_NI * _NI[_qp];
  _D_GB[_qp] = D_GB;
  _h[_qp] = ShamNeedlemann::h_psi(psi);
  _b_sat[_qp] = (1. / std::sqrt(M_PI * Nmax_NI * _NI[_qp]));
  _E_GB[_qp] = E_GB;
  _G_GB[_qp] = G_GB;
  _eta_sliding[_qp] = eta_sliding;
  _sigma_0[_qp] = sigma_0;
  _thickness[_qp] = thickness;
  if (_scale_stiffness)
    _thickness[_qp] *= _JxW_ref / _JxW[_qp];
  _beta_exponent[_qp] = beta_exponent;
  _n_exponent[_qp] = n_exponent;
}

void
GBCavitation::getInitPropertyValuesFromParams(Real & FN_NI,
                                              Real & Nmax_NI,
                                              Real & a0,
                                              Real & b0,
                                              Real & psi,
                                              Real & D_GB,
                                              Real & E_GB,
                                              Real & G_GB,
                                              Real & eta_sliding,
                                              Real & sigma_0,
                                              Real & thickness,
                                              Real & beta_exponent,
                                              Real & n_exponent) const
{
  FN_NI = getParam<Real>("FN_NI");
  Nmax_NI = getParam<Real>("Nmax_NI");
  a0 = getParam<Real>("a0");
  b0 = getParam<Real>("b0");
  psi = getParam<Real>("psi_degree");
  D_GB = getParam<Real>("D_GB");
  E_GB = getParam<Real>("E_GB");
  G_GB = getParam<Real>("G_GB");
  eta_sliding = getParam<Real>("eta_sliding");
  sigma_0 = getParam<Real>("sigma_0");
  thickness = getParam<Real>("interface_thickness");
  beta_exponent = getParam<Real>("beta_exponent");
  n_exponent = getParam<Real>("n_exponent");
}

void
GBCavitation::getInitPropertyValuesFromUO(Real & FN_NI,
                                          Real & Nmax_NI,
                                          Real & a0,
                                          Real & b0,
                                          Real & psi,
                                          Real & D_GB,
                                          Real & E_GB,
                                          Real & G_GB,
                                          Real & eta_sliding,
                                          Real & sigma_0,
                                          Real & thickness,
                                          Real & beta_exponent,
                                          Real & n_exponent) const
{
  const std::map<std::string, Real> prop_map =
      _GBCavitationBoundaryPropertyUO->getPropertyMap(_current_elem->id(), _current_side);

  auto ptr = prop_map.find("FN_NI");
  if (ptr != prop_map.end())
    FN_NI = ptr->second;
  else
    mooseError("can't find FN_NI ");

  ptr = prop_map.find("Nmax_NI");
  if (ptr != prop_map.end())
    Nmax_NI = ptr->second;
  else
    mooseError("can't find Nmax_NI ");

  ptr = prop_map.find("a0");
  if (ptr != prop_map.end())
    a0 = ptr->second;
  else
    mooseError("can't find a0 ");

  ptr = prop_map.find("b0");
  if (ptr != prop_map.end())
    b0 = ptr->second;
  else
    mooseError("can't find b0 ");

  ptr = prop_map.find("psi_degree");
  if (ptr != prop_map.end())
    psi = ptr->second;
  else
    mooseError("can't find psi_degree ");

  ptr = prop_map.find("D_GB");
  if (ptr != prop_map.end())
    D_GB = ptr->second;
  else
    mooseError("can't find D_GB ");

  ptr = prop_map.find("E_GB");
  if (ptr != prop_map.end())
    E_GB = ptr->second;
  else
    mooseError("can't find E_GB ");

  ptr = prop_map.find("G_GB");
  if (ptr != prop_map.end())
    G_GB = ptr->second;
  else
    mooseError("can't find G_GB ");

  ptr = prop_map.find("eta_sliding");
  if (ptr != prop_map.end())
    eta_sliding = ptr->second;
  else
    mooseError("can't find eta_sliding ");

  ptr = prop_map.find("sigma_0");
  if (ptr != prop_map.end())
    sigma_0 = ptr->second;
  else
    mooseError("can't find sigma_0 ");

  ptr = prop_map.find("thickness");
  if (ptr != prop_map.end())
    thickness = ptr->second;
  else
    mooseError("can't find thickness ");

  ptr = prop_map.find("n_exponent");
  if (ptr != prop_map.end())
    n_exponent = ptr->second;
  else
    mooseError("can't find n_exponent ");

  ptr = prop_map.find("beta_exponent");
  if (ptr != prop_map.end())
    beta_exponent = ptr->second;
  else
    mooseError("can't find beta_exponent ");
}
