#pragma once

#include "MathUtils.h"

#include "v1dot.h"

v1dot::v1dot(const std::string &fname, const qFun &q, const Real &D,
             const Real &theta_time_integration)
    : nlFunRateNoHistory(fname, theta_time_integration), _q(q), _D(D) {}

Real v1dot::computeValue(const io_maps_type &x, const io_maps_type &params,
                         const io_maps_type &x_old) const {
  UNUSED(x_old);
  Real Tn = getValueFromMap(x, "Tn", "v1dot::x");
  Real q = _q.computeValue(x, params, x_old);

  Real v = 8 * _pi * _D * Tn / q;
  if (std::isfinite(v))
    return v;
  else
    mooseError(_f_name + " value is not finite " + std::to_string(v) + ": " +
               std::to_string(q));
  return 0;
}

Real v1dot::dV1dqFun(const Real &q, const Real &Tn) const {
  return -8 * _pi * _D * Tn / (q * q);
}

v1dot::io_maps_type v1dot::computeVarGradient(const io_maps_type &x,
                                              const io_maps_type &params,
                                              const io_maps_type &x_old) const {
  UNUSED(x_old);
  Real Tn = getValueFromMap(x, "Tn", "v1dot::x");

  Real f = 8 * _pi * _D * Tn;
  io_maps_type df_dx = {{"Tn", 8 * _pi * _D}};
  Real g = _q.computeValue(x, params, x_old);
  io_maps_type dg_dx = _q.computeVarGradient(x, params, x_old);

  return f_divided_g_D(f, df_dx, g, dg_dx);
}

v1dot::io_maps_type
v1dot::computeParamGradient(const io_maps_type &x, const io_maps_type &params,
                            const io_maps_type &x_old) const {
  UNUSED(x_old);
  Real Tn = getValueFromMap(x, "Tn", "v1dot::x");

  Real f = 8 * _pi * _D * Tn;
  io_maps_type df_dp = _empty_map;
  Real g = _q.computeValue(x, params, x_old);
  io_maps_type dg_dp = _q.computeVarGradient(x, params, x_old);

  return f_divided_g_D(f, df_dp, g, dg_dp);
}
