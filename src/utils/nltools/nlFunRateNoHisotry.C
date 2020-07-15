#pragma once

#include "nlFunRateNoHistory.h"

nlFunRateNoHistory::nlFunRateNoHistory(const std::string &fname,
                                       const Real &theta_time_integration)
    : nlFunBase(fname, requiredVar(), requiredParam()),
      _theta_implicit(1 - theta_time_integration),
      _theta_explicit(theta_time_integration) {}

Real nlFunRateNoHistory::computeTimeIntegral(const io_maps_type &x,
                                             const io_maps_type &params,
                                             const io_maps_type &x_old,
                                             const Real &dt) const {

  Real V_dot_impl = 0;
  Real V_dot_expl = 0;
  computeValue(x, params, x_old);

  if (std::abs(_theta_implicit) > 1e-6)
    V_dot_impl = _theta_implicit * computeValue(x, params, x_old);
  if (std::abs(_theta_explicit) > 1e-6)
    V_dot_expl = _theta_explicit * computeValue(x_old, params, x_old);

  return (V_dot_impl + V_dot_expl) * dt;
}

nlFunRateNoHistory::io_maps_type
nlFunRateNoHistory::computeTimeIntegralVarGradient(const io_maps_type &x,
                                                   const io_maps_type &params,
                                                   const io_maps_type &x_old,
                                                   const Real &dt) const {

  io_maps_type dVdot_dx = {};
  if (std::abs(_theta_implicit) > 1e-6) {
    dVdot_dx = computeVarGradient(x, params, x_old);
    dVdot_dx = chain(dt * _theta_implicit, dVdot_dx);
  }

  return dVdot_dx;
}

nlFunRateNoHistory::io_maps_type
nlFunRateNoHistory::computeTimeIntegralParamGradient(const io_maps_type &x,
                                                     const io_maps_type &params,
                                                     const io_maps_type &x_old,
                                                     const Real &dt) const {

  io_maps_type dVdot_dp_impl = {};
  io_maps_type dVdot_dp_expl = {};
  if (std::abs(_theta_implicit) > 1e-6) {
    dVdot_dp_impl = computeParamGradient(x, params, x_old);
    dVdot_dp_impl = chain(dt * _theta_implicit, dVdot_dp_impl);
  }

  if (std::abs(_theta_explicit) > 1e-6) {
    dVdot_dp_expl = computeParamGradient(x_old, params, x_old);
    dVdot_dp_expl = chain(dt * _theta_explicit, dVdot_dp_expl);
  }

  return sumD(dVdot_dp_impl, dVdot_dp_expl);
}
