#pragma once

#include "nlFunRate.h"

nlFunRate::nlFunRate(const std::string &fname, const std::string &var_name,
                     const Real &theta_time_integration)
    : nlFunRateNoHistory(fname, theta_time_integration), _var_name(var_name) {}

Real nlFunRate::computeTimeIntegral(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old,
                                    const Real &dt) const {

  Real V_old =
      getValueFromMap(x_old, _var_name, "nlFunRate::computeTimeIntegral.x_old");

  Real V_int =
      V_old + nlFunRateNoHistory::computeTimeIntegral(x, params, x_old, dt);
  return V_int;
}

// nlFunRate::io_maps_type nlFunRate::computeTimeIntegralVarGradient(
//     const io_maps_type &x, const io_maps_type &params,
//     const io_maps_type &x_old, const Real &dt) const {
//
//   io_maps_type dVdot_dx = computeVarGradient(x, params, x_old);
//   dVdot_dx = chain(dt, dVdot_dx);
//
//   return dVdot_dx;
// }
//
// nlFunRate::io_maps_type nlFunRate::computeTimeIntegralParamGradient(
//     const io_maps_type &x, const io_maps_type &params,
//     const io_maps_type &x_old, const Real &dt) const {
//
//   io_maps_type dVdot_dp = computeParamGradient(x, params, x_old);
//   dVdot_dp = chain(dt, dVdot_dp);
//
//   return dVdot_dp;
// }
