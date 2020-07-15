//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBCavitation_smooth.h"
#include "libmesh/quadrature.h"

registerMooseObject("DeerApp", GBCavitation_smooth);

InputParameters GBCavitation_smooth::validParams() {
  InputParameters params = CZMMaterialBasePK1::validParams();
  params.addParam<Real>("a0", 5e-5, "intial caivity half radius");
  params.addParam<Real>("b0", 6e-2, "intial caivity half spacing");
  params.addParam<Real>("Nmax_NI", 1000, "saturation cavity half spacing");
  params.addParam<Real>("FN_NI", 2e4, "normalized nucleation rate constant");
  params.addParam<Real>("S_thr", 0,
                        "the threshould value at which nucleation begins");
  params.addParam<Real>("E_interface", 150e3, "interface Young modulus");
  params.addParam<Real>("E_penalty", 10, "innerpenetration peanlty");
  params.addParam<Real>("G_interface", 52.63e3, "interface shear modulus");
  params.addParam<Real>("beta_exponent", 2, "nucleation stress exponent");
  params.addParam<Real>("D_gb", 1e-15, "gran boundary diffusion coefficient");
  params.addParam<Real>("n_exponent", 5, "power law creep exponent");
  params.addParam<Real>("psi_angle", 75,
                        "equilibrium cavity tip half-angle [degree]");
  params.addParam<Real>("sigma_0", 200, "traction normalization parameter");
  params.addParam<Real>("D_thr", 0.95,
                        "the maximum damage before the elemet is killed");
  params.addParam<Real>("max_allowed_opening_traction", 10000,
                        "After reaching this traction value the element is "
                        "considered unstable, and therefore is killed");
  params.addParam<Real>("max_allowed_damage_rate", 1e-4,
                        "After reaching this damage rate the element is "
                        "considered unstable, and therefore is killed");

  params.addParam<Real>("min_allowed_residual_stiffness", 1e-6,
                        "the minimum allowed stiffness after element kill");
  params.addParam<Real>("traction_mean_decay_time_factor", 0.5,
                        "the maximum damage bifore the elemet is killed");
  params.addParam<Real>("eta_sliding", 1e6, "sliding viscosity");
  params.addParam<bool>("use_old_avg_prop", true,
                        "if true (default) old matrial property will "
                        "be used for the average interface quantities");
  params.addParam<int>("linesearch_type", 1, "0 no linesearch, 1 backtracking");
  params.addParam<bool>("snes_use_FD", false, "use Finte difference jacobian");
  params.addParam<Real>("snes_abs_tol", 1e-4,
                        "abs tolerandce for material point snes ");
  params.addParam<Real>("snes_rel_tol", 1e-6,
                        "realtive tolerandce for material point snes ");
  params.addParam<Real>("snes_step_tol", 1e-6,
                        "step tolerance for convergence ");
  params.addParam<Real>("snes_max_iteration", 50,
                        "step tolerance for convergence ");
  params.addParam<int>(
      "vdot_type", 3,
      "1 VL, 2 VH,  3(default)  max(VL,VH), 4 average (VL+VH)/2");
  params.addParam<int>(
      "triaxial_vdot_active", 0,
      "false deactivate the V1H and/or V2H, 0 (default) means inactive");
  params.addParam<int>("vdot_max_type", 0,
                       "0(default) hard max, 1 softMax, 2 softMaxMin, 3 choose "
                       "on previous step vdot max");
  params.addParam<Real>("vdot_smooth_factor", 20,
                        " soft max smoothing paramter");
  params.addParam<bool>("use_substep", true, "allow to use substep");
  params.addParam<bool>("force_substep", false, "allow to use substep");
  params.addParam<bool>("use_LM", true,
                        "use lagrangian multipliares to impse constriant on "
                        "the values of a and b");
  params.addParam<bool>("cavity_nucleation_on", true, "tunrs on bbot equation");
  params.addParam<bool>("cavity_growth_on", true, "turn on adot equation");
  params.addParam<bool>("triaxial_cavity_growth_on", true,
                        "enable triaxial cavity growth");
  params.addParam<Real>("min_allowed_residual_life", 50,
                        "enable triaxial cavity growth");
  params.addParam<Real>(
      "theta_time_integration", 0.5,
      "theta level of explicitness int eh theta method: 1 forward euler (fully "
      "ecplicit), 0 backward euler (fully implicit order 1), "
      "0.5 crank-nicolson  (fully implicit order 2), all other values generate "
      "order 1 methods with different degrees of explicintess");
  params.addParam<bool>("give_up_qp", false,
                        "true if you want to use previous stae solution if a "
                        "qp does not ocnverge");
  params.addParam<unsigned int>("max_substep_cuts", 5,
                                "the maximum number of substep cuts");
  return params;
}

GBCavitation_smooth::GBCavitation_smooth(const InputParameters &parameters)
    : CZMMaterialBasePK1(parameters),
      _GBCavitationBoundaryPropertyUO(
          parameters.isParamSetByUser("GBCavitationBoundaryPropertyUO")
              ? &getUserObjectByName<GBCavitationBoundaryPropertyUO>(
                    getParam<UserObjectName>("GBCavitationBoundaryPropertyUO"))
              : nullptr),
      _displacement_jump_dot(
          declareProperty<RealVectorValue>("displacement_jump_dot")),
      _displacement_jump_dot_old(
          getMaterialPropertyOld<RealVectorValue>("displacement_jump_dot")),
      _E_penalty(getParam<Real>("E_penalty")),
      _beta_exponent(getParam<Real>("beta_exponent")),
      _n_exponent(getParam<Real>("n_exponent")),
      _alpha_n(3. / (_n_exponent * 2.)), _a0(declareProperty<Real>("a0")),
      _a0_old(getMaterialPropertyOld<Real>("a0")),
      _b0(declareProperty<Real>("b0")),
      _b0_old(getMaterialPropertyOld<Real>("b0")),
      _NI(declareProperty<Real>("NI")),
      _NI_old(getMaterialPropertyOld<Real>("NI")),
      _FN(declareProperty<Real>("FN")),
      _FN_old(getMaterialPropertyOld<Real>("FN")),
      _D_gb(declareProperty<Real>("D_gb")),
      _D_gb_old(getMaterialPropertyOld<Real>("D_gb")),
      _b_saturation(declareProperty<Real>("b_saturation")),
      _b_saturation_old(getMaterialPropertyOld<Real>("b_saturation")),
      _E_interface(declareProperty<Real>("E_interface")),
      _E_interface_old(getMaterialPropertyOld<Real>("E_interface")),
      _G_interface(declareProperty<Real>("G_interface")),
      _G_interface_old(getMaterialPropertyOld<Real>("G_interface")),
      _eta_sliding(declareProperty<Real>("eta_sliding")),
      _eta_sliding_old(getMaterialPropertyOld<Real>("_eta_sliding")),
      _interface_thickness(declareProperty<Real>("interface_thickness")),
      _interface_thickness_old(
          getMaterialPropertyOld<Real>("interface_thickness")),
      _sigma_0(declareProperty<Real>("sigma_0")),
      _sigma_0_old(getMaterialPropertyOld<Real>("sigma_0")),
      _S_thr(declareProperty<Real>("S_thr")),
      _S_thr_old(getMaterialPropertyOld<Real>("S_thr")),
      _h(declareProperty<Real>("h")), _h_old(getMaterialPropertyOld<Real>("h")),
      _residualScaleFactors("residaul_scale_factor"),
      /* interface properties*/
      _a(declareProperty<Real>("average_cavity_radii")),
      _a_old(getMaterialPropertyOld<Real>("average_cavity_radii")),
      _b(declareProperty<Real>("average_cavity_spacing")),
      _b_old(getMaterialPropertyOld<Real>("average_cavity_spacing")),
      _D(declareProperty<Real>("interface_damage")),
      _D_old(getMaterialPropertyOld<Real>("interface_damage")),
      _D_dot(declareProperty<Real>("interface_damage_rate")),
      _D_dot_old(getMaterialPropertyOld<Real>("interface_damage_rate")),
      _residual_life(declareProperty<Real>("residual_life")),
      _residual_life_old(getMaterialPropertyOld<Real>("residual_life")),
      _VL1_dot(declareProperty<Real>("VL1_dot")),
      _VL2_dot(declareProperty<Real>("VL2_dot")),
      _VH1_dot(declareProperty<Real>("VH1_dot")),
      _VH2_dot(declareProperty<Real>("VH2_dot")),
      _Vdot(declareProperty<Real>("Vdot")),
      _VLdot(declareProperty<Real>("VLdot")),
      _VLdot_old(getMaterialPropertyOld<Real>("VLdot")),
      _VHdot(declareProperty<Real>("VHdot")),
      _VHdot_old(getMaterialPropertyOld<Real>("VHdot")),
      _nucleation_above_threshold(
          declareProperty<bool>("nucleation_above_threshold")),
      _nucleation_above_threshold_old(
          getMaterialPropertyOld<bool>("nucleation_above_threshold")),
      /* bulk material related properties */
      _use_old_avg_prop(getParam<bool>("use_old_avg_prop")),
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
      _avg_mises_stress(declareProperty<Real>("avg_sVM_interface")),
      _avg_mises_stress_old(getMaterialPropertyOld<Real>("avg_sVM_interface")),

      _avg_hyd_stress(declareProperty<Real>("avg_sH_interface")),
      _avg_hyd_stress_old(getMaterialPropertyOld<Real>("avg_sH_interface")),

      _avg_eq_strain_rate(declareProperty<Real>("avg_eceq_interface")),
      _avg_eq_strain_rate_old(
          getMaterialPropertyOld<Real>("avg_eceq_interface")),

      _accumulated_eq_strain(declareProperty<Real>("accumulated_eq_strain")),
      _accumulated_eq_strain_old(
          getMaterialPropertyOld<Real>("accumulated_eq_strain")),
      _interface_triaxiality(declareProperty<Real>("interface_triaxiality")),
      /* failure related properties */
      _elem_failed(declareProperty<bool>("elem_failed")),
      _elem_failed_old(getMaterialPropertyOld<bool>("elem_failed")),
      _D_thr(getParam<Real>("D_thr")),
      _max_allowed_opening_traction(
          getParam<Real>("max_allowed_opening_traction")),
      _max_allowed_damage_rate(getParam<Real>("max_allowed_damage_rate")),
      _min_allowed_residual_life(getParam<Real>("min_allowed_residual_life")),
      _traction_mean_decay_time_factor(
          getParam<Real>("traction_mean_decay_time_factor")),
      _min_allowed_residual_stiffness(
          getParam<Real>("min_allowed_residual_stiffness")),
      _time_at_failure(declareProperty<Real>("time_at_failure")),
      _time_at_failure_old(getMaterialPropertyOld<Real>("time_at_failure")),
      _traction_at_failure(
          declareProperty<RealVectorValue>("traction_at_failure")),
      _traction_at_failure_old(
          getMaterialPropertyOld<RealVectorValue>("traction_at_failure")),

      // _decay_exhausted(declareProperty<bool>("decay_exhausted")),
      // _decay_exhausted_old(getMaterialPropertyOld<bool>("decay_exhausted")),
      _du_at_failure(declareProperty<RealVectorValue>("du_at_failure")),
      _du_at_failure_old(
          getMaterialPropertyOld<RealVectorValue>("du_at_failure")),
      _K_at_failure(declareProperty<RealVectorValue>("stiffness_at_failure")),
      _K_at_failure_old(
          getMaterialPropertyOld<RealVectorValue>("stiffness_at_failure")),
      _mysnes_abs_tol(getParam<Real>("snes_abs_tol")),
      _mysnes_rel_tol(getParam<Real>("snes_rel_tol")),
      _mysnes_step_tol(getParam<Real>("snes_step_tol")),
      _mysnes_max_iteration(getParam<Real>("snes_max_iteration")),
      _use_substep(getParam<bool>("use_substep")),
      _max_substep_cuts(getParam<unsigned int>("max_substep_cuts")),
      _force_substep(getParam<bool>("force_substep")),
      _use_LM(getParam<bool>("use_LM")), _n_equation(_use_LM ? 6 : 3),
      _give_up_qp(getParam<bool>("give_up_qp")) {
  // initialize petsc required variables
  _q_ierr = SNESCreate(PETSC_COMM_SELF, &_q_snes);
  SNESSetOptionsPrefix(_q_snes, "GBCavitation_smooth");
  _q_ierr = VecCreate(PETSC_COMM_SELF, &_q_x);

  _q_ierr = VecSetSizes(_q_x, PETSC_DECIDE, _n_equation);
  _q_ierr = VecSetFromOptions(_q_x);
  _q_ierr = VecDuplicate(_q_x, &_q_r);

  /* Create Jacobian matrix data structure */
  _q_ierr = MatCreate(PETSC_COMM_SELF, &_q_J);
  _q_ierr =
      MatSetSizes(_q_J, PETSC_DECIDE, PETSC_DECIDE, _n_equation, _n_equation);
  _q_ierr = MatSetFromOptions(_q_J);
  _q_ierr = MatSetUp(_q_J);

  /*  Set function evaluation routine and vector */
  SNESSetFunction(_q_snes, _q_r, FormFunction1, &_ctx);

  /* Set Jacobian matrix data structure and Jacobian evaluation routine */
  if (!parameters.get<bool>("snes_use_FD"))
    SNESSetJacobian(_q_snes, _q_J, _q_J, FormJacobian1, &_ctx);
  else
    SNESSetJacobian(_q_snes, _q_J, _q_J, SNESComputeJacobianDefault, &_ctx);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */
  if (parameters.get<int>("linesearch_type") == 0)
    setNewton(/*linesearch_on=*/false);
  if (parameters.get<int>("linesearch_type") == 1)
    setNewton(/*linesearch_on=*/true);
}

GBCavitation_smooth::~GBCavitation_smooth() { _q_ierr = SNESDestroy(&_q_snes); }

void GBCavitation_smooth::initQpStatefulProperties() {
  CZMMaterialBasePK1::initQpStatefulProperties();

  // init real vector value properties
  for (unsigned int i = 0; i < 3; i++) {
    _traction_at_failure[_qp](i) = 0;
    _du_at_failure[_qp](i) = 0;
    _K_at_failure[_qp](i) = 0;
    _displacement_jump_dot[_qp](i) = 0;
  }

  _elem_failed[_qp] = false;
  _time_at_failure[_qp] = 0;

  _VHdot[_qp] = 0;
  _VLdot[_qp] = 0;
  _nucleation_above_threshold[_qp] = false;
  _avg_mises_stress[_qp] = 0;
  _avg_hyd_stress[_qp] = 0;
  _avg_eq_strain_rate[_qp] = 0;
  _accumulated_eq_strain[_qp] = 0;
  _residual_life[_qp] = 1e6;

  InitGBCavitationParamsAndProperties();
}

void GBCavitation_smooth::InitGBCavitationParamsAndProperties() {
  Real FN_NI, Nmax_NI, a0, b0, psi, D_gb, E_interface, G_interface, eta_sliding,
      sigma_0, S_thr;
  if (_GBCavitationBoundaryPropertyUO == nullptr)
    getInitPropertyValuesFromParams(FN_NI, Nmax_NI, a0, b0, psi, D_gb,
                                    E_interface, G_interface, eta_sliding,
                                    sigma_0, S_thr);
  else
    getInitPropertyValuesFromUO(FN_NI, Nmax_NI, a0, b0, psi, D_gb, E_interface,
                                G_interface, eta_sliding, sigma_0, S_thr);

  psi *= libMesh::pi / 180.;

  _a0[_qp] = a0;
  _b0[_qp] = b0;
  _a[_qp] = a0;
  _b[_qp] = b0;
  _D[_qp] = a0 / b0;
  _D_dot[_qp] = 0;
  _NI[_qp] = (1. / (libMesh::pi * b0 * b0));
  _FN[_qp] = FN_NI * _NI[_qp];
  _D_gb[_qp] = D_gb;
  _h[_qp] = (1. / (1. + std::cos(psi)) - std::cos(psi) / 2.);
  _b_saturation[_qp] = b0 * std::pow(Nmax_NI, -0.5);
  _E_interface[_qp] = E_interface;
  _G_interface[_qp] = G_interface;
  _eta_sliding[_qp] = eta_sliding;
  _interface_thickness[_qp] = _b_saturation[_qp] * 6;
  _sigma_0[_qp] = sigma_0;
  _S_thr[_qp] = S_thr;
}

void GBCavitation_smooth::getInitPropertyValuesFromParams(
    Real &FN_NI, Real &Nmax_NI, Real &a0, Real &b0, Real &psi, Real &D_gb,
    Real &E_interface, Real &G_interface, Real &eta_sliding, Real &sigma_0,
    Real &S_thr) const {
  a0 = getParam<Real>("a0");
  b0 = getParam<Real>("b0");
  FN_NI = getParam<Real>("FN_NI");
  Nmax_NI = getParam<Real>("Nmax_NI");
  a0 = getParam<Real>("a0");
  b0 = getParam<Real>("b0");
  psi = getParam<Real>("psi_angle");
  D_gb = getParam<Real>("D_gb");
  E_interface = getParam<Real>("E_interface");
  G_interface = getParam<Real>("G_interface");
  eta_sliding = getParam<Real>("eta_sliding");
  sigma_0 = getParam<Real>("sigma_0");
  S_thr = _pars.isParamSetByUser("S_thr") ? getParam<Real>("S_thr")
                                          : 1. / getParam<Real>("FN_NI");
}

void GBCavitation_smooth::getInitPropertyValuesFromUO(
    Real &FN_NI, Real &Nmax_NI, Real &a0, Real &b0, Real &psi, Real &D_gb,
    Real &E_interface, Real &G_interface, Real &eta_sliding, Real &sigma_0,
    Real &S_thr) const {
  const std::map<std::string, Real> prop_map =
      _GBCavitationBoundaryPropertyUO->getPropertyMap(_current_elem->id(),
                                                      _current_side);

  auto ptr = prop_map.find("Nmax_NI");
  if (ptr != prop_map.end())
    Nmax_NI = ptr->second;
  else
    mooseError("can't find Nmax_NI ");

  ptr = prop_map.find("FN_NI");
  if (ptr != prop_map.end())
    FN_NI = ptr->second;
  else
    mooseError("can't find FN_NI ");

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

  ptr = prop_map.find("D_gb");
  if (ptr != prop_map.end())
    D_gb = ptr->second;
  else
    mooseError("can't find D_gb ");

  ptr = prop_map.find("psi_angle");
  if (ptr != prop_map.end())
    psi = ptr->second;
  else
    mooseError("can't find psi_angle ");

  ptr = prop_map.find("E_interface");
  if (ptr != prop_map.end())
    E_interface = ptr->second;
  else
    mooseError("can't find E_interface ");

  ptr = prop_map.find("G_interface");
  if (ptr != prop_map.end())
    G_interface = ptr->second;
  else
    mooseError("can't find G_interface ");

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

  ptr = prop_map.find("S_thr");
  if (ptr != prop_map.end())
    S_thr = ptr->second;
  else
    mooseError("can't find S_thr ");
}

void GBCavitation_smooth::computeTractionIncrementAndDerivatives() {
  _displacement_jump_dot[_qp] = _displacement_jump_inc[_qp] / _dt;
  if (_t_step > 0) {
    if (!_elem_failed_old[_qp]) {
      computeAverageBulkPorperties();
      _a0[_qp] = _a0_old[_qp];
      _b0[_qp] = _b0_old[_qp];
      _a[_qp] = _a_old[_qp];
      _b[_qp] = _b_old[_qp];
      _D[_qp] = _D_old[_qp];
      _D_dot[_qp] = _D_dot_old[_qp];
      _D_gb[_qp] = _D_gb_old[_qp];
      _h[_qp] = _h_old[_qp];
      _E_interface[_qp] = _E_interface_old[_qp];
      _G_interface[_qp] = _G_interface_old[_qp];
      _eta_sliding[_qp] = _eta_sliding_old[_qp];
      _interface_thickness[_qp] = _interface_thickness_old[_qp];
      _b_saturation[_qp] = _b_saturation_old[_qp];
      _sigma_0[_qp] = _sigma_0_old[_qp];
      _S_thr[_qp] = _S_thr_old[_qp];
      _FN[_qp] = _FN_old[_qp];
      _NI[_qp] = _NI_old[_qp];

      GBCavitationNLSystem GBNLsystem(
          "GBNLSYS", _D_gb[_qp], _h[_qp], _n_exponent, _E_interface[_qp],
          _E_penalty, _G_interface[_qp], _eta_sliding[_qp],
          _interface_thickness[_qp], _beta_exponent, _a0_old[_qp], _b0_old[_qp],
          _b_saturation[_qp], _sigma_0[_qp], _S_thr[_qp], _FN[_qp],
          getParam<int>("vdot_max_type"), getParam<int>("vdot_type"),
          getParam<int>("triaxial_vdot_active"),
          getParam<Real>("vdot_smooth_factor"),
          getParam<bool>("cavity_nucleation_on"),
          getParam<bool>("cavity_growth_on"),
          getParam<bool>("triaxial_cavity_growth_on"),
          getParam<Real>("theta_time_integration"));

      /* init old property map */
      nlFunBase::io_maps_type x_sol_real;
      nlFunBase::io_maps_type x_old = xoldFromOld();

      /* init paramters  map */
      nlFunBase::io_maps_type NL_params = {
          {"eqedotc", _avg_eq_strain_rate_old[_qp]},
          {"eqec", _accumulated_eq_strain_old[_qp]},
          {"sVM", _avg_mises_stress_old[_qp]},
          {"sH", _avg_hyd_stress_old[_qp]},
          {"un", _displacement_jump[_qp](0)},
          {"us1", _displacement_jump[_qp](1)},
          {"us2", _displacement_jump[_qp](2)},
          {"un_dot", _displacement_jump_dot[_qp](0)},
          {"un_dot_old", _displacement_jump_dot_old[_qp](0)},
          {"us1_dot", _displacement_jump_dot[_qp](1)},
          {"us2_dot", _displacement_jump_dot[_qp](2)},
          {"un_old", _displacement_jump_old[_qp](0)},
          {"a0", _a0_old[_qp]},
      };

      _nucleation_above_threshold[_qp] = _nucleation_above_threshold_old[_qp];
      if (!_nucleation_above_threshold_old[_qp]) {
        _nucleation_above_threshold[_qp] =
            GBNLsystem.nucleationAboveThreshold(NL_params, x_old);
      }
      if (!_nucleation_above_threshold[_qp])
        NL_params["nucleation_above_threshold"] = 0;
      else {
        NL_params["nucleation_above_threshold"] = 1;
      }

      _interface_triaxiality[_qp] = 0;
      if (_avg_mises_stress[_qp] != 0)
        _interface_triaxiality[_qp] =
            _avg_hyd_stress[_qp] / _avg_mises_stress[_qp];
      /* initialize petsc context */
      _ctx.dt = _dt;
      _ctx.GBNLsysstem_pt = &GBNLsystem;
      _ctx.my_xold = &x_old;
      _ctx.my_params = &NL_params;
      _residual_scale_factors =
          _residualScaleFactors.computeVarScaleFactor(x_old);
      _ctx.scale_factor_pt = &_residual_scale_factors;
      _ctx.n_eq_pt = &_n_equation;

      /* initialize petsc intial guess */
      initSnesGuess(x_old, x_old, GBNLsystem);

      /* PETSC solve */
      _q_ierr = SNESSolve(_q_snes, NULL, _q_x);

      /* check convergence */
      bool not_converged = true;
      if (!_force_substep) {
        bool not_converged = !checkCavitationConvergence(
            NL_params, x_old, _dt, x_sol_real, GBNLsystem);

        if (!not_converged) {
          update_Dtn_dUN(NL_params, x_old, _dt, GBNLsystem);
          GBNLsystem.returnVolumeRate(x_sol_real, NL_params, x_old,
                                      _VL1_dot[_qp], _VL2_dot[_qp],
                                      _VH1_dot[_qp], _VH2_dot[_qp], _Vdot[_qp],
                                      _VLdot[_qp], _VHdot[_qp]);
          updatedStateVarFromRealSolution(x_sol_real, _dt);
        }
      }

      bool fail_while_substep = false;
      if (_force_substep || (not_converged && _use_substep))
        not_converged = !substepFun(x_sol_real, fail_while_substep, GBNLsystem);

      if (not_converged)
        throw MooseException("GB_cavitation failed!");

      // if (_give_up_qp) {
      //   Moose::out << " giving up intehration of element "
      //              << _current_elem->id() << " qp " << _qp << "\n";
      //   _throw_exception_mp[_qp] = false;
      //   copyOldSolution();
      // } else {
      //   // setNewton(/*linesearch = */ false);
      //   // writeDataToFile();
      //   return;
      // }

      /* if we are here the material point solution has converged*/
      // _throw_exception_mp[_qp] = false;
      // return;
      // if (!fail_while_substep) {
      //
      //   /* solution achieved, material property have already been updated in
      //   the
      //    * check*/
      //   _throw_exception_mp[_qp] = false;
      //   // decoupeldShearTraction(_dt);
      //
      // } else { /* failed_faile_substep*/
      //          /* traction and derivtives already updated */
      //   // mooseWarning("fail while substep");
      // }

    } // end element not failed
    else if (_elem_failed_old[_qp]) {
      tractionDeacy(_time_at_failure_old[_qp], _a_old[_qp], _b_old[_qp],
                    _D_old[_qp], _D_dot_old[_qp], _residual_life_old[_qp],
                    _traction_at_failure_old[_qp], _du_at_failure_old[_qp],
                    _K_at_failure_old[_qp]);
      _VLdot[_qp] = 0;
      _VHdot[_qp] = 0;

      _VL1_dot[_qp] = 0;
      _VL2_dot[_qp] = 0;
      _VH1_dot[_qp] = 0;
      _VH2_dot[_qp] = 0;
      _Vdot[_qp] = 0;

    }      /*endelemtn failed */
  } else { /*_t_step == 0*/
    InitGBCavitationParamsAndProperties();

    _D_dot[_qp] = 0;
    _VLdot[_qp] = 0;
    _VHdot[_qp] = 0;

    _VL1_dot[_qp] = 0;
    _VL2_dot[_qp] = 0;
    _VH1_dot[_qp] = 0;
    _VH2_dot[_qp] = 0;
    _Vdot[_qp] = 0;

    _time_at_failure[_qp] = 0;
    _elem_failed[_qp] = false;
    // _decay_exhausted[_qp] = false;
    for (unsigned int i = 0; i < 3; i++) {
      _traction[_qp](i) = 0;
      _traction_at_failure[_qp](i) = 0;
      _du_at_failure[_qp](i) = 0;
      _K_at_failure[_qp](i) = 0;
    }
    _avg_mises_stress[_qp] = 0;
    _avg_hyd_stress[_qp] = 0;
    _avg_eq_strain_rate[_qp] = 0;
    _residual_life[_qp] = 1e6;
  }
  _traction_inc[_qp] = _traction[_qp] - _traction_old[_qp];
}

/* ------------------------------------------------------------------- */
/*
   FormFunction1 - Evaluates nonlinear function, F(x).

   Input Parameters:
.  snes - the SNES context
.  x    - input vector
.  _ctx  - optional user-defined context

   Output Parameter:
.  f - function vector
 */
PetscErrorCode GBCavitation_smooth::FormFunction1(SNES q_snes, Vec x, Vec f,
                                                  void *ctx) {
  ApplicationCtx *user = (ApplicationCtx *)ctx;
  nlFunBase::io_maps_type *xold_petsc = user->my_xold;
  nlFunBase::io_maps_type *params_petsc = user->my_params;
  Real dt_pestc = user->dt;
  GBCavitationNLSystem *GBNLsystem_petsc = user->GBNLsysstem_pt;
  const nlFunBase::io_maps_type *sF = user->scale_factor_pt;
  // const int *n_eq = user->n_eq_pt;
  const int n_var = 3;
  PetscErrorCode ierr;
  const PetscScalar *xxl;
  PetscScalar *ff;

  /*
   Get pointers to vector data.
      - For default PETSc vectors, VecGetArray() returns a pointer to
        the data array.  Otherwise, the routine is implementation dependent.
      - You MUST call VecRestoreArray() when you no longer need access to
        the array.
   */
  nlFunBase::io_maps_type curr_x;
  VecGetArrayRead(x, &xxl);
  curr_x["a"] = xxl[0];
  curr_x["b"] = xxl[1];
  curr_x["Tn"] = xxl[2];
  // curr_x["Ts1"] = xxl[3];
  // curr_x["Ts2"] = xxl[4];
  VecGetArray(f, &ff);
  std::vector<Real> S(3, 0);
  std::vector<std::string> vname{"a", "b", "Tn", "Ts1", "Ts2"};
  for (unsigned int i = 0; i < 3; i++)
    S[i] = sF->find(vname[i])->second;

  DenseVector<Real> x_comp = GBNLsystem_petsc->computeSystemValue(
      curr_x, *params_petsc, *xold_petsc, dt_pestc);
  /* Compute function */
  for (unsigned int i = 0; i < 3; i++)
    ff[i] = (xxl[i] - x_comp(i)) / S[i];

  // if (n_var > 3) {
  //   unsigned int idx_a = 0;
  //   unsigned int idx_b = 1;
  //   unsigned int idx_l1 = 3;
  //   unsigned int idx_l2 = 4;
  //   unsigned int idx_l3 = 5;
  //   /* constraint b>a -> g = a-b */
  //   /* min(lambda, -g)*/
  //   Real l1 = xxl[idx_l1];
  //   Real g1 = xxl[0] - xxl[1];
  //   ff[idx_l1] = l1 > -g1 ? -g1 : l1;
  //
  //   // /* constraint a>_a0 */
  //   Real l2 = xxl[idx_l2];
  //   Real g2 = params_petsc->find("a0")->second - xxl[0];
  //   ff[idx_l2] = l2 > -g2 ? -g2 : l2;
  //
  //   // /* constraint a>_a0 */
  //   Real l3 = xxl[idx_l3];
  //   Real g3 = xxl[1] - xold_petsc->find("b")->second;
  //   ff[idx_l3] = l3 > -g3 ? -g3 : l3;
  //
  //   ff[0] = ff[0] + xxl[idx_l1] - xxl[idx_l2];
  //   ff[1] = ff[1] - xxl[idx_l1] + xxl[idx_l3];
  //
  //   // ff[0] = ff[0] + xxl[5] - xxl[7];
  //   // ff[1] = ff[1] - xxl[5] - xxl[6];
  // }
  /* Restore vectors */
  VecRestoreArrayRead(x, &xxl);
  VecRestoreArray(f, &ff);
  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormJacobian1 - Evaluates Jacobian matrix.

   Input Parameters:
.  q_snes - the SNES context
.  x - input vector
.  ctx - bringing in required data

   Output Parameters:
.  jac - Jacobian matrix
.  B - optionally different preconditioning matrix
.  flag - flag indicating matrix structure
*/
PetscErrorCode GBCavitation_smooth::FormJacobian1(SNES q_snes, Vec x, Mat jac,
                                                  Mat B, void *ctx) {
  ApplicationCtx *user = (ApplicationCtx *)ctx;
  nlFunBase::io_maps_type *xold_petsc = user->my_xold;
  nlFunBase::io_maps_type *params_petsc = user->my_params;
  Real dt_petsc = user->dt;
  GBCavitationNLSystem *GBNLsystem_petsc = user->GBNLsysstem_pt;
  const nlFunBase::io_maps_type *sF_petsc = user->scale_factor_pt;
  // const int *n_eq = user->n_eq_pt;
  PetscInt n_var = 3;
  PetscErrorCode ierr;
  const PetscScalar *xxl;
  PetscScalar A[n_var * n_var];
  PetscInt idx[n_var];
  for (int i = 0; i < n_var; i++)
    idx[i] = i;

  /*
     Get pointer to vector data, create map for NL equation and then restore
  */
  nlFunBase::io_maps_type curr_x;
  VecGetArrayRead(x, &xxl);
  curr_x["a"] = xxl[0];
  curr_x["b"] = xxl[1];
  curr_x["Tn"] = xxl[2];
  // curr_x["Ts1"] = xxl[3];
  // curr_x["Ts2"] = xxl[4];

  std::vector<Real> S(3, 0);
  std::vector<std::string> vname{"a", "b", "Tn", "Ts1", "Ts2"};
  for (unsigned int i = 0; i < 3; i++)
    S[i] = sF_petsc->find(vname[i])->second;

  /*
     Compute Jacobian entries and insert into matrix.
  */

  DenseMatrix<Real> J = GBNLsystem_petsc->computeSystemVarJacobian(
      curr_x, *params_petsc, *xold_petsc, dt_petsc, /*wrt_xNL=*/true);

  for (unsigned int i = 0; i < n_var * n_var; i++)
    A[i] = 0;

  for (unsigned int eq = 0; eq < 3; eq++)
    for (unsigned int var = 0; var < 3; var++) {
      Real deq_dvar = -J(eq, var);
      if (eq == var)
        deq_dvar += 1;
      A[eq * n_var + var] = deq_dvar / S[eq];
    }

  // if (n_var > 3) {
  //   unsigned int idx_a = 0;
  //   unsigned int idx_b = 1;
  //   unsigned int idx_l1 = 3;
  //   unsigned int idx_l2 = 4;
  //   unsigned int idx_l3 = 5;
  //
  //   Real l1 = xxl[idx_l1];
  //   Real g1 = std::abs(xxl[0]) - std::abs(xxl[1]);
  //
  //   Real l2 = xxl[idx_l2];
  //   Real g2 = params_petsc->find("a0")->second - xxl[0];
  //
  //   Real l3 = xxl[idx_l3];
  //   Real g3 = xxl[1] - xold_petsc->find("b")->second;
  //
  //   A[n_var * idx_a + idx_l1] = 1;
  //   A[n_var * idx_a + idx_l2] = -1.;
  //
  //   A[n_var * idx_b + idx_l1] = -1;
  //   A[n_var * idx_b + idx_l3] = 1;
  //
  //   A[n_var * idx_l1 + idx_a] = l1 > -g1 ? -1 : 0;
  //   A[n_var * idx_l1 + idx_b] = l1 > -g1 ? 1 : 0;
  //   A[n_var * idx_l1 + idx_l1] = l1 > -g1 ? 0 : 1;
  //
  //   A[n_var * idx_l2 + idx_a] = l2 > -g2 ? 1 : 0;
  //   A[n_var * idx_l2 + idx_l2] = l2 > -g2 ? 0 : 1;
  //
  //   A[n_var * idx_l3 + idx_b] = l3 > -g3 ? -1 : 0;
  //   A[n_var * idx_l3 + idx_l3] = l3 > -g3 ? 0 : 1;
  // }

  VecRestoreArrayRead(x, &xxl);

  MatSetValues(B, n_var, idx, n_var, idx, A, INSERT_VALUES);

  /*
     Assemble matrix
  */
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  if (jac != B) {
    MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
  }
  return 0;
}

/* mostly for debug, write data to text file */
void GBCavitation_smooth::writeDataToFile() const {
  std::string fname = "elem_" + std::to_string(_current_elem->id()) + ".data";
  std::ofstream outfile;
  outfile.open(fname, std::ios_base::app);
  /*
    format
    "t:"   "_qp "   " _dt: "   " T_old: "   " a_old "   " b_old "   " u_old: "
    " u_dot: "   " sVM "   " sH "   " edoteq "   " eeq "
    */
  outfile << _t - _dt << " " << _qp << " " << _dt << " "
          << _traction_old[_qp](0) << " " << _traction_old[_qp](1) << " "
          << _traction_old[_qp](2) << " " << _a_old[_qp] << " " << _b_old[_qp]
          << " " << _displacement_jump[_qp](0) << " "
          << _displacement_jump[_qp](1) << " " << _displacement_jump[_qp](2)
          << " " << _displacement_jump_old[_qp](0) << " "
          << _displacement_jump_old[_qp](1) << " "
          << _displacement_jump_old[_qp](2) << " "
          << _displacement_jump_dot[_qp](0) << " "
          << _displacement_jump_dot[_qp](1) << " "
          << _displacement_jump_dot[_qp](2) << " " << _avg_mises_stress[_qp]
          << " " << _avg_hyd_stress[_qp] << " " << _avg_eq_strain_rate[_qp]
          << " " << _accumulated_eq_strain[_qp] << " " << _t_step << std::endl;
}

void GBCavitation_smooth::setNewton(const bool &linesearch_on) {
  /* set solver type*/
  _q_ierr = SNESGetType(_q_snes, &_q_snestype);
  // _q_ierr = SNESSetType(_q_snes, SNESNEWTONLS);
  _q_ierr = SNESSetType(_q_snes, SNESNEWTONLS);
  /* set solver tolerances*/

  _q_ierr = SNESSetTolerances(_q_snes, _mysnes_abs_tol, _mysnes_rel_tol,
                              _mysnes_step_tol, _mysnes_max_iteration, -1);
  /* set line search */
  _q_ierr = SNESGetLineSearch(_q_snes, &_q_linesearch);

  if (!linesearch_on)
    _q_ierr = SNESLineSearchSetType(_q_linesearch, SNESLINESEARCHBASIC);
  else
    _q_ierr = SNESLineSearchSetType(_q_linesearch, SNESLINESEARCHBT);

  /* KSP and PC */
  _q_ierr = SNESGetKSP(_q_snes, &_q_ksp);
  _q_ierr = KSPGetPC(_q_ksp, &_q_pc);
  _q_ierr = PCSetType(_q_pc, PCLU);

  // KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,PetscReal
  // dtol,PetscInt maxits)
  _q_ierr =
      KSPSetTolerances(_q_ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 20);

  /* apply changes to snes solver */
  _q_ierr = SNESSetFromOptions(_q_snes);
}

bool GBCavitation_smooth::substepFun(nlFunBase::io_maps_type &x_sol_real,
                                     bool &fail_while_substep,
                                     GBCavitationNLSystem &GBNLsystem) {
  Real substep_dt = _dt / 2.;
  Real substep_dt_old = substep_dt;
  unsigned int n_time_cut = 1;
  Real time_last_solution = 0;

  nlFunBase::io_maps_type x_init = xoldFromOld();
  x_init["VL_dot"] = _VLdot_old[_qp];
  x_init["VH_dot"] = _VHdot_old[_qp];

  nlFunBase::io_maps_type x_old = x_init;
  nlFunBase::io_maps_type x_older;

  x_sol_real = x_old;
  nlFunBase::io_maps_type NL_params, params_old;
  fail_while_substep = false;
  // _nucleation_above_threshold[_qp] =

  do {
    NL_params = prepareParamsSubstep(time_last_solution, substep_dt);
    if (!_nucleation_above_threshold[_qp])
      NL_params["nucleation_above_threshold"] = 0;
    else
      NL_params["nucleation_above_threshold"] = 1;

    prepareSolverContext(substep_dt, x_old, NL_params, GBNLsystem);
    initSnesGuess(x_sol_real, x_old, GBNLsystem);
    _q_ierr = SNESSolve(_q_snes, NULL, _q_x);

    if (checkCavitationConvergence(NL_params, x_old, substep_dt, x_sol_real,
                                   GBNLsystem)) {
      /* set x_old to new solution*/
      GBNLsystem.returnVolumeRate(x_sol_real, NL_params, x_old, _VL1_dot[_qp],
                                  _VL2_dot[_qp], _VH1_dot[_qp], _VH2_dot[_qp],
                                  _Vdot[_qp], _VLdot[_qp], _VHdot[_qp]);

      /*update time*/
      time_last_solution += substep_dt;
      /*check if failure is predicted*/
      fail_while_substep = false;
      if ((x_sol_real.find("a")->second / x_sol_real.find("b")->second) >=
          _D_thr)
        fail_while_substep = true;

      if (!fail_while_substep) {
        if (time_last_solution != _dt) {
          /*prepare for next step*/
          /*xsol -> xold*/
          x_old = x_sol_real;
          x_old["VL_dot"] = _VLdot_old[_qp];
          x_old["VH_dot"] = _VHdot_old[_qp];

          /* check substep size*/
          if (time_last_solution + substep_dt >= _dt)
            substep_dt = _dt - time_last_solution;
        }
      }
    } else {
      /* reduce time */
      substep_dt = substep_dt / 2.;
      n_time_cut += 1;
      /* reset guess to previous known solution*/
      x_sol_real = x_old;
      x_sol_real.erase("VL_dot");
      x_sol_real.erase("VH_dot");
    }
  } while (n_time_cut < _max_substep_cuts && time_last_solution != _dt &&
           !fail_while_substep);
  if (time_last_solution == _dt) {
    /* SUBSTEP STANDARD CONVERGENCE*/
    update_Dtn_dUN(NL_params, x_old, substep_dt, GBNLsystem);
    updatedStateVarFromRealSolution(x_sol_real, _dt);
    return true;
  } else if (fail_while_substep) {
    /*FAIL WHILE SUBSTEPPING, NEED TO UPDATE STATE VAR AND CALL TRACTCTION
     * DECAY*/
    updatedStateVarFromRealSolution(x_sol_real, substep_dt);
    GBNLsystem.returnVolumeRate(x_sol_real, NL_params, x_old, _VL1_dot[_qp],
                                _VL2_dot[_qp], _VH1_dot[_qp], _VH2_dot[_qp],
                                _Vdot[_qp], _VLdot[_qp], _VHdot[_qp]);

    _elem_failed[_qp] = true;

    tractionDeacy(_t - (_dt - time_last_solution), _a[_qp], _b[_qp], _D[_qp],
                  _D_dot[_qp], _residual_life[_qp], _traction_at_failure[_qp],
                  _du_at_failure[_qp], _K_at_failure[_qp]);

    return true;
  } else {
    /*SUBSTEP DID NOT CONVERGE*/
    return false;
  }
}

void GBCavitation_smooth::prepareSolverContext(
    const Real &substep_dt, nlFunBase::io_maps_type &x_old,
    nlFunBase::io_maps_type &NL_params, GBCavitationNLSystem &GBNLsystem) {
  _ctx.dt = substep_dt;
  _ctx.GBNLsysstem_pt = &GBNLsystem;
  _ctx.my_xold = &x_old;
  _ctx.my_params = &NL_params;
  _residual_scale_factors = _residualScaleFactors.computeVarScaleFactor(x_old);
  _ctx.scale_factor_pt = &_residual_scale_factors;
  _ctx.n_eq_pt = &_n_equation;
}

nlFunBase::io_maps_type
GBCavitation_smooth::prepareParamsSubstep(const Real &time_last_solution,
                                          const Real &current_dt) {
  nlFunBase::io_maps_type NL_params;
  Real total_dt = time_last_solution + current_dt;

  NL_params = {
      {"eqedotc", _avg_eq_strain_rate_old[_qp]},
      {"eqec", _accumulated_eq_strain_old[_qp]},
      {"sVM", _avg_mises_stress_old[_qp]},
      {"sH", _avg_hyd_stress_old[_qp]},
      {"un", _displacement_jump_old[_qp](0) +
                 _displacement_jump_dot[_qp](0) * total_dt},
      {"us1", _displacement_jump_old[_qp](1) +
                  _displacement_jump_dot[_qp](1) * total_dt},
      {"us2", _displacement_jump_old[_qp](2) +
                  _displacement_jump_dot[_qp](2) * total_dt},
      {"un_dot", _displacement_jump_dot[_qp](0)},
      {"un_dot_old", _displacement_jump_dot_old[_qp](0)},
      {"us1_dot", _displacement_jump_dot[_qp](1)},
      {"us2_dot", _displacement_jump_dot[_qp](2)},
      {"un_old", _displacement_jump_old[_qp](0) +
                     _displacement_jump_dot[_qp](0) * time_last_solution},
      {"a0", _a0_old[_qp]},
  };

  return NL_params;
}

void GBCavitation_smooth::initSnesGuess(
    const nlFunBase::io_maps_type &x_Real,
    const nlFunBase::io_maps_type &x_Real_old,
    const GBCavitationNLSystem &GBNLsystem) {
  PetscScalar *xx;
  nlFunBase::io_maps_type xNL =
      GBNLsystem.getNLValueFromRealValue(x_Real, x_Real_old);
  VecGetArray(_q_x, &xx);
  xx[0] = xNL.find("a")->second;
  xx[1] = xNL.find("b")->second;
  xx[2] = xNL.find("Tn")->second;
  // xx[3] = xNL.find("Ts1")->second;
  // xx[4] = xNL.find("Ts2")->second;
  if (_use_LM)
    for (unsigned int i = 3; i < _n_equation; i++)
      xx[i] = 0.;

  VecRestoreArray(_q_x, &xx);
}

void GBCavitation_smooth::updateLastSolution(
    nlFunBase::io_maps_type &x_last_solution) const {
  const PetscScalar *xxc;

  VecGetArrayRead(_q_x, &xxc);
  x_last_solution["a"] = _a[_qp];
  x_last_solution["b"] = _b[_qp];
  x_last_solution["Tn"] = _traction[_qp](0);
  // x_last_solution["Ts1"] = _traction[_qp](1);
  // x_last_solution["Ts2"] = _traction[_qp](2);
  VecRestoreArrayRead(_q_x, &xxc);
}

bool GBCavitation_smooth::checkCavitationConvergence(
    const nlFunBase::io_maps_type &NL_params,
    const nlFunBase::io_maps_type &x_old, const Real &dt_local,
    nlFunBase::io_maps_type &x_sol_real,
    const GBCavitationNLSystem &GBNLsystem) {
  SNESGetConvergedReason(_q_snes, &_q_reason);
  if (_q_reason < 0 || _q_reason > 4) {
    // PetscPrintf(PETSC_COMM_SELF, "%s: \n",
    // SNESConvergedReasons[_q_reason]);
    return false;
  }

  x_sol_real =
      getRealSolutionFromNLSolution(NL_params, x_old, dt_local, GBNLsystem);
  Real a = x_sol_real.find("a")->second;
  Real b = x_sol_real.find("b")->second;

  if (b < 0 || a > b)
    return false;

  if (b > _b_old[_qp]) {
    mooseWarning("current b grater than b_old: b= " + std::to_string(b) +
                 " b_old= " + std::to_string(_b_old[_qp]) +
                 " diff= " + std::to_string(std::abs(_b_old[_qp] - b)));
    return false;
  }

  return true;
}

nlFunBase::io_maps_type GBCavitation_smooth::getRealSolutionFromNLSolution(
    const nlFunBase::io_maps_type &NL_params,
    const nlFunBase::io_maps_type &x_old, const Real &dt_local,
    const GBCavitationNLSystem &GBNLsystem) const {
  nlFunBase::io_maps_type xNL = copyNLsolutionToMap();
  nlFunBase::io_maps_type x_Real =
      GBNLsystem.getRealValueFromNLValue(xNL, x_old);

  return GBNLsystem.computeSystemRealValue(x_Real, NL_params, x_old, dt_local);
}

nlFunBase::io_maps_type GBCavitation_smooth::copyNLsolutionToMap() const {
  const PetscScalar *xxc;
  nlFunBase::io_maps_type my_NL_solution = {};
  VecGetArrayRead(_q_x, &xxc);
  my_NL_solution["a"] = xxc[0];
  my_NL_solution["b"] = xxc[1];
  my_NL_solution["Tn"] = xxc[2];
  // my_NL_solution["Ts1"] = xxc[3];
  // my_NL_solution["Ts2"] = xxc[4];
  VecRestoreArrayRead(_q_x, &xxc);

  return my_NL_solution;
}

void GBCavitation_smooth::updatedStateVarFromRealSolution(
    const nlFunBase::io_maps_type &x_Real, const Real &dt_curr) {
  _a[_qp] = x_Real.find("a")->second;
  _b[_qp] = x_Real.find("b")->second;
  _traction[_qp](0) = x_Real.find("Tn")->second;
  decoupeldShearTraction(dt_curr);

  _D[_qp] = _a[_qp] / _b[_qp];
  _elem_failed[_qp] = false;

  _D_dot[_qp] = (_D[_qp] - _D_old[_qp]) / dt_curr;
  _residual_life[_qp] = 1e6;
  if (_D_dot[_qp] > 0)
    _residual_life[_qp] = (1. - _D[_qp]) / _D_dot[_qp];

  Real traction_rate = (_traction[_qp](0) - _traction_old[_qp](0)) / _dt;

  if (_D[_qp] >= _D_thr) /*" exceeded maximum Damage, marked as failed"*/
    _elem_failed[_qp] = true;

  else if ((_traction[_qp](0)) >= _max_allowed_opening_traction &&
           _t_step > 10) /*" exceeded maximum Traction, marked as failed"*/
    _elem_failed[_qp] = true;

  else if (_traction[_qp](0) > 0 && _D_dot[_qp] > 0 &&
           _residual_life[_qp] <=
               _min_allowed_residual_life) /*" below minimum reidual life,
                                               marked as failed"*/
    _elem_failed[_qp] = true;

  /*the following helps to prevent bad compression cases causing oscillations*/
  else if (((_residual_life[_qp] < _min_allowed_residual_life) ||
            (traction_rate * _min_allowed_residual_life + _traction[_qp](0) <=
             -_max_allowed_opening_traction)) &&
           (_t_step > 10) && (_traction[_qp](0) < 0) &&
           (_displacement_jump_dot[_qp](0) > 0))
    _elem_failed[_qp] = true;

  else
    _elem_failed[_qp] = false;

  for (unsigned int i = 0; i < 3; i++) {
    _traction_at_failure[_qp](i) = _traction[_qp](i);
    _du_at_failure[_qp](i) = _displacement_jump_old[_qp](i) +
                             _displacement_jump_dot[_qp](i) * dt_curr;
    _K_at_failure[_qp](i) =
        std::abs(_traction_at_failure[_qp](i) / _du_at_failure[_qp](i));
  }
  _time_at_failure[_qp] = _t - (_dt - dt_curr);
  if (_elem_failed[_qp])
    _D[_qp] = _D_thr;
}

nlFunBase::io_maps_type GBCavitation_smooth::xoldFromOld() const {
  nlFunBase::io_maps_type x_old = {{"a", _a_old[_qp]},
                                   {"b", _b_old[_qp]},
                                   {"Tn", _traction_old[_qp](0)},
                                   {"Ts1", _traction_old[_qp](1)},
                                   {"Ts2", _traction_old[_qp](2)},
                                   {"VL_dot", _VLdot_old[_qp]},
                                   {"VH_dot", _VHdot_old[_qp]}};
  return x_old;
}

void GBCavitation_smooth::decoupeldShearTraction(const Real &dt) {
  Real a = _a_old[_qp];
  Real b = _b_old[_qp];
  Real S = _interface_thickness[_qp] / (_G_interface[_qp] * (1. - (a / b)));
  Real TS1, TS2, dTS1_duS1, dTS2_duS2;
  Real eta = _eta_sliding[_qp];
  if (a / b > 0.5)
    eta = 2 * _eta_sliding[_qp] * (1. - a / b);
  _traction[_qp](1) =
      eta * _displacement_jump_dot[_qp](1) +
      std::exp(-dt / (eta * S)) *
          (_traction_old[_qp](1) - eta * _displacement_jump_dot[_qp](1));

  _dtraction_djump[_qp](1, 0) = 0;

  _dtraction_djump[_qp](1, 1) = eta / dt - std::exp(-dt / (eta * S)) * eta / dt;
  _dtraction_djump[_qp](1, 2) = 0;

  _traction[_qp](2) =
      eta * _displacement_jump_dot[_qp](2) +
      std::exp(-dt / (eta * S)) *
          (_traction_old[_qp](2) - eta * _displacement_jump_dot[_qp](2));

  _dtraction_djump[_qp](2, 0) = 0;
  _dtraction_djump[_qp](2, 1) = 0;
  _dtraction_djump[_qp](2, 2) = eta / dt - std::exp(-dt / (eta * S)) * eta / dt;
}

void GBCavitation_smooth::computeAverageBulkPorperties() {
  // compute average Von Mises Stress;
  RankTwoTensor dev_stress = _stress_master[_qp].deviatoric();
  _avg_mises_stress[_qp] =
      std::sqrt(1.5 * dev_stress.doubleContraction(dev_stress));
  dev_stress = _stress_slave[_qp].deviatoric();
  _avg_mises_stress[_qp] +=
      std::sqrt(1.5 * dev_stress.doubleContraction(dev_stress));
  _avg_mises_stress[_qp] /= 2.;

  // compute average hydrostatic stress;
  _avg_hyd_stress[_qp] =
      (_stress_master[_qp].trace() + _stress_slave[_qp].trace()) / 6.;

  // compute equivalent inelastic strain rate;
  RankTwoTensor strain_rate =
      (_inelastic_strain_master[_qp] - _inelastic_strain_master_old[_qp]) / _dt;
  _avg_eq_strain_rate[_qp] =
      std::sqrt(2.0 / 3.0 * strain_rate.doubleContraction(strain_rate));
  strain_rate =
      (_inelastic_strain_slave[_qp] - _inelastic_strain_slave_old[_qp]) / _dt;
  _avg_eq_strain_rate[_qp] +=
      std::sqrt(2.0 / 3.0 * strain_rate.doubleContraction(strain_rate));
  _avg_eq_strain_rate[_qp] /= 2.;

  _accumulated_eq_strain[_qp] =
      _accumulated_eq_strain_old[_qp] + _avg_eq_strain_rate[_qp] * _dt;
}
