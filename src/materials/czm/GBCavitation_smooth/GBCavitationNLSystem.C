#pragma once

#include "GBCavitationNLSystem.h"

GBCavitationNLSystem::GBCavitationNLSystem(
    const std::string &fname, const Real &D_gb, const Real &h,
    const Real &n_exponent, const Real &E_interface, const Real &E_penalty,
    const Real &G_interface, const Real &eta_sliding,
    const Real &interface_thickness, const Real &beta_exponent, const Real &a0,
    const Real &b0, const Real &b_saturation, const Real &sigma_0,
    const Real &S_thr, const Real &FN, const int &vdot_max_type,
    const int &vdot_type, const int &triaxial_vdot_active,
    const Real &vdot_smooth_factor, const bool &cavity_nucleation_on,
    const bool &cavity_growth_on, const bool &triaxial_cavity_growth_on,
    const Real &theta_time_integration)
    : nlFunBase(fname),
      _eqTests("eqTests", a0, b0, sigma_0, S_thr, beta_exponent, b_saturation,
               E_interface, cavity_nucleation_on, cavity_growth_on,
               triaxial_cavity_growth_on),
      _fab("fab"), _faL("faL", D_gb), _fLow("fLow", _faL, _fab),
      _qLow("qLow", _fLow), _qHigh("qHigh", _fab),
      _v1L("v1Low", _qLow, D_gb, theta_time_integration),
      _v1H("v1High", _qHigh, D_gb, theta_time_integration),
      _v2L("v2Low", n_exponent, h, theta_time_integration, _eqTests),
      _v2H("v2High", n_exponent, h, theta_time_integration, _eqTests),
      _vdot("vDot", _v1L, _v1H, _v2L, _v2H, _eqTests, vdot_max_type, vdot_type,
            triaxial_vdot_active, vdot_smooth_factor, theta_time_integration),
      _Ts1dot("Ts1dot", "us1", "Ts1", G_interface, interface_thickness,
              eta_sliding, theta_time_integration),
      _Ts2dot("Ts2dot", "us2", "Ts2", G_interface, interface_thickness,
              eta_sliding, theta_time_integration),
      _TNdot("TNdot", "un", "Tn", E_interface, interface_thickness, _vdot,
             E_penalty, theta_time_integration),
      _adot("adot", _vdot, h, a0, _eqTests, theta_time_integration),
      _bdot("bdot", beta_exponent, b0, b_saturation, S_thr, sigma_0, FN,
            _eqTests, theta_time_integration),
      _aComp("aComp", _adot, a0, theta_time_integration),
      _bComp("bComp", _bdot, theta_time_integration),
      _prepSNVar("prepVars", a0), _vdot_type(vdot_type),
      _triaxial_vdot_active(triaxial_vdot_active) {}

DenseVector<Real> GBCavitationNLSystem::computeSystemValue(
    const io_maps_type &x_NL_in, const io_maps_type &params,
    const io_maps_type &x_old, const Real &dt) const {
  io_maps_type my_x = _prepSNVar.computeNewX(x_NL_in, x_old);

  DenseVector<Real> computed_values(3);

  io_maps_type computed_x = computeSystemRealValue(my_x, params, x_old, dt);

  io_maps_type x_NL_out = _prepSNVar.XNLFromComputedX(computed_x, x_old);

  computed_values(0) = getValueFromMap(x_NL_out, "a");
  computed_values(1) = getValueFromMap(x_NL_out, "b");
  computed_values(2) = getValueFromMap(x_NL_out, "Tn");
  // computed_values(3) = getValueFromMap(x_NL_out, "Ts1");
  // computed_values(4) = getValueFromMap(x_NL_out, "Ts2");

  return computed_values;
}

GBCavitationNLSystem::io_maps_type GBCavitationNLSystem::computeSystemRealValue(
    const io_maps_type &x_Real, const io_maps_type &params,
    const io_maps_type &x_old, const Real &dt) const {

  io_maps_type computed_x = _empty_map;
  computed_x["a"] = _aComp.computeTimeIntegral(x_Real, params, x_old, dt);
  computed_x["b"] = _bComp.computeTimeIntegral(x_Real, params, x_old, dt);
  computed_x["Tn"] = _TNdot.computeTimeIntegral(x_Real, params, x_old, dt);
  // computed_x["Ts1"] = _Ts1dot.computeTimeIntegral(x_Real, params, x_old, dt);
  // computed_x["Ts2"] = _Ts2dot.computeTimeIntegral(x_Real, params, x_old, dt);

  return computed_x;
}

void GBCavitationNLSystem::varGradientHelper(const io_maps_type &var_grad,
                                             const io_maps_type &myx_dx,
                                             DenseMatrix<Real> &J,
                                             const unsigned int &idx) const {
  const std::vector<std::string> v_name{"a", "b", "Tn", "Ts1", "Ts2"};
  for (unsigned int i = 0; i < 3; i++) {
    J(idx, i) = 0;
    if (checkMapKeyExist(var_grad, v_name[i]))
      J(idx, i) = getValueFromMap(var_grad, v_name[i]) *
                  getValueFromMap(myx_dx, v_name[i]);
  }
}

DenseMatrix<Real> GBCavitationNLSystem::computeSystemVarJacobian(
    const io_maps_type &xNL, const io_maps_type &params,
    const io_maps_type &x_old, const Real &dt, const bool &wrt_XNL) const {
  DenseMatrix<Real> computed_var_gradient(3, 3);
  io_maps_type var_grad;
  io_maps_type my_x = _prepSNVar.computeNewX(xNL, x_old);
  io_maps_type myx_dx = my_x;
  if (wrt_XNL)
    myx_dx = _prepSNVar.computeNewXVarGradient(xNL, x_old);
  else
    for (auto &it : myx_dx)
      it.second = 1;

  var_grad = _aComp.computeTimeIntegralVarGradient(my_x, params, x_old, dt);
  varGradientHelper(var_grad, myx_dx, computed_var_gradient, 0);

  var_grad = _bComp.computeTimeIntegralVarGradient(my_x, params, x_old, dt);
  varGradientHelper(var_grad, myx_dx, computed_var_gradient, 1);

  var_grad = _TNdot.computeTimeIntegralVarGradient(my_x, params, x_old, dt);
  varGradientHelper(var_grad, myx_dx, computed_var_gradient, 2);

  // var_grad = _Ts1dot.computeTimeIntegralVarGradient(my_x, params, x_old, dt);
  // varGradientHelper(var_grad, myx_dx, computed_var_gradient, 3);
  //
  // var_grad = _Ts2dot.computeTimeIntegralVarGradient(my_x, params, x_old, dt);
  // varGradientHelper(var_grad, myx_dx, computed_var_gradient, 4);
  return computed_var_gradient;
}

std::map<std::string, DenseVector<Real>>
GBCavitationNLSystem::computeSystemParamGradient(const io_maps_type &xNL,
                                                 const io_maps_type &params,
                                                 const io_maps_type &x_old,
                                                 const Real &dt) const {
  io_maps_type my_x = _prepSNVar.computeNewX(xNL, x_old);
  // io_maps_type myx_dx = x;
  // if (wrt_XNL)
  //   myx_dx = _prepSNVar.computeNewXVarGradient(xNL, x_old);
  // else
  //   for (auto &it : myx_dx)
  //     it->second = 1;

  std::map<std::string, DenseVector<Real>> dx_dparam;

  io_maps_type da_dp =
      _aComp.computeTimeIntegralParamGradient(my_x, params, x_old, dt);
  io_maps_type db_dp =
      _bComp.computeTimeIntegralParamGradient(my_x, params, x_old, dt);
  io_maps_type dTN_dp =
      _TNdot.computeTimeIntegralParamGradient(my_x, params, x_old, dt);
  // io_maps_type dTS1_dp =
  //     _Ts1dot.computeTimeIntegralParamGradient(my_x, params, x_old, dt);
  // io_maps_type dTS2_dp =
  //     _Ts2dot.computeTimeIntegralParamGradient(my_x, params, x_old, dt);

  for (const auto &p : params) {
    DenseVector<Real> dx_dp(3, 0);

    if (checkMapKeyExist(da_dp, p.first))
      dx_dp(0) =
          getValueFromMap(da_dp, p.first); // * getValueFromMap(myx_dx, "a");

    if (checkMapKeyExist(db_dp, p.first))
      dx_dp(1) =
          getValueFromMap(db_dp, p.first); // * getValueFromMap(myx_dx, "b");

    if (checkMapKeyExist(dTN_dp, p.first))
      dx_dp(2) =
          getValueFromMap(dTN_dp, p.first); // * getValueFromMap(myx_dx, "Tn");

    // if (checkMapKeyExist(dTS1_dp, p.first))
    //   dx_dp(3) = getValueFromMap(dTS1_dp,
    //                              p.first); // * getValueFromMap(myx_dx,
    //                              "Ts1");
    //
    // if (checkMapKeyExist(dTS2_dp, p.first))
    //   dx_dp(4) = getValueFromMap(dTS2_dp,
    //                              p.first); // * getValueFromMap(myx_dx,
    //                              "Ts2");

    dx_dparam[p.first] = dx_dp;
  }

  return dx_dparam;
}

nlFunBase::io_maps_type
GBCavitationNLSystem::getRealValueFromNLValue(const io_maps_type &x_NL,
                                              const io_maps_type &x_old) const {
  return _prepSNVar.computeNewX(x_NL, x_old);
}

nlFunBase::io_maps_type
GBCavitationNLSystem::getNLValueFromRealValue(const io_maps_type &xcomp,
                                              const io_maps_type &x_old) const {
  return _prepSNVar.XNLFromComputedX(xcomp, x_old);
}

void GBCavitationNLSystem::returnVolumeRate(const io_maps_type &x,
                                            const io_maps_type &params,
                                            const io_maps_type &x_old,
                                            Real &VL1, Real &VL2, Real &VH1,
                                            Real &VH2, Real &Vdot, Real &VLdot,
                                            Real &VHdot) {

  VL1 = 0;
  VL2 = 0;
  VH1 = 0;
  VH2 = 0;
  if (_vdot_type == 1 || _vdot_type == 3) {
    VL1 = _v1L.computeValue(x, params, x_old);
    if (_eqTests.triaxialGrowthIsActive(x, params, x_old)) {
      VL2 = _v2L.computeValue(x, params, x_old);
    }
  }
  VLdot = VL1 + VL2;
  if (_vdot_type == 2 || _vdot_type == 3) {
    VH1 = _v1H.computeValue(x, params, x_old);
    if (_eqTests.triaxialGrowthIsActive(x, params, x_old)) {
      VH2 = _v2H.computeValue(x, params, x_old);
    }
  }
  VHdot = VH1 + VH2;
  Vdot = _vdot.computeValue(x, params, x_old);
}
