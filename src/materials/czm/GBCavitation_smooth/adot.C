#pragma once

#include "adot.h"

adot::adot(const std::string &fname, const vdot &vdotFun, const Real &h,
           const Real &a0, const equationTests &eqTest,
           const Real &theta_time_integration)
    : nlFunRate(fname, "a", theta_time_integration), _vdot(vdotFun), _h(h),
      _a0(a0), _eqTest(eqTest) {}

Real adot::computeValue(const io_maps_type &x, const io_maps_type &params,
                        const io_maps_type &x_old) const {

  UNUSED(x_old);

  if (!_eqTest.growthIsActive(x, params, x_old))
    return 0;

  Real a = getValueFromMap(x, "a", "adot::computeValue.x");
  Real vdot_value = _vdot.computeValue(x, params, x_old);
  Real a_dot = vdot_value / (4. * _pi * _h * std::pow(a, 2.));

  return a_dot;
}

adot::io_maps_type adot::computeVarGradient(const io_maps_type &x,
                                            const io_maps_type &params,
                                            const io_maps_type &x_old) const {
  UNUSED(x_old);
  if (!_eqTest.growthIsActive(x, params, x_old))
    return _empty_map;

  Real f = _vdot.computeValue(x, params, x_old);
  io_maps_type df_dx = _vdot.computeVarGradient(x, params, x_old);

  Real a = getValueFromMap(x, "a", "adot::computeValue.x");
  Real g = (4. * _pi * _h * std::pow(a, 2.));
  io_maps_type dg_dx = {{"a", 8. * _pi * _h * a}};

  return f_divided_g_D(f, df_dx, g, dg_dx);
}

adot::io_maps_type adot::computeParamGradient(const io_maps_type &x,
                                              const io_maps_type &params,
                                              const io_maps_type &x_old) const {
  return _empty_map;
}

Real adot::dfdaFun(const io_maps_type &x, const io_maps_type &params,
                   const io_maps_type &x_old) const {
  Real a = getValueFromMap(x, "a", "adot::dfdaFun::x");
  Real dfda =
      -_vdot.computeValue(x, params, x_old) / (2. * _pi * _h * std::pow(a, 3.));
  return dfda;
}

Real adot::dfdvdotFun(const io_maps_type &x) const {
  Real a = getValueFromMap(x, "a", "adot::dfdaFun::x");
  Real dfdvdot = 1. / (4. * _pi * _h * std::pow(a, 2.));
  return dfdvdot;
}
