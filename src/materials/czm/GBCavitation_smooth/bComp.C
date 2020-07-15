#pragma once

#include "bComp.h"

bComp::bComp(const std::string &fname, const bdot &bFun,
             const Real &theta_time_integration)
    : nlFunRate(fname, "b", theta_time_integration), _bdot(bFun),
      _soft_min("soft_min", 100, true, false, true) {}

Real bComp::computeTimeIntegral(const io_maps_type &x,
                                const io_maps_type &params,
                                const io_maps_type &x_old,
                                const Real &dt) const {

  // io_maps_type b_values = {
  //     {"bmax", b0.computeValue(x, params, x_old)},
  //     {"bComp", _bdot.computeTimeIntegral(x, params, x_old, dt)}};
  //
  // Real f = _soft_min.computeValue(b_values, params, x_old);

  // nlFunConstant b0("bmax", getValueFromMap(x_old, "b"));
  // Real f = std::min(b0.computeValue(x, params, x_old),
  //                   std::abs(_bdot.computeTimeIntegral(x, params, x_old,
  //                   dt)));

  Real f = _bdot.computeTimeIntegral(x, params, x_old, dt);

  if (std::isfinite(f))
    return f;
  else
    mooseError(_f_name + "value is not finite " + std::to_string(f));
  return 0;
}

bComp::io_maps_type bComp::computeTimeIntegralVarGradient(
    const io_maps_type &x, const io_maps_type &params,
    const io_maps_type &x_old, const Real &dt) const {

  // nlFunConstant b0("bmax", getValueFromMap(x_old, "b"));

  // Real Bint = _bdot.computeTimeIntegral(x, params, x_old, dt);
  //
  // // if (std::abs(Bint) < b0.computeValue(x, params, x_old))
  // return chain(std::copysign(1., Bint),
  //              _bdot.computeTimeIntegralVarGradient(x, params, x_old, dt));
  // else
  //   return b0.computeVarGradient(x, params, x_old);
  // io_maps_type b_values = {
  //     {"bmax", b0.computeValue(x, params, x_old)},
  //     {"bComp", _bdot.computeTimeIntegral(x, params, x_old, dt)}};
  //
  // io_maps_type b0_grad = b0.computeVarGradient(x, params, x_old);
  // io_maps_type b_grad =
  //     _bdot.computeTimeIntegralVarGradient(x, params, x_old, dt);
  // map_of_io_maps_type b_values_grad = {{"bmax", &b0_grad}, {"bComp",
  // &b_grad}};

  // return _soft_min.computeVarGradient(b_values, b_values_grad, params);

  return _bdot.computeTimeIntegralVarGradient(x, params, x_old, dt);
}

bComp::io_maps_type bComp::computeTimeIntegralParamGradient(
    const io_maps_type &x, const io_maps_type &params,
    const io_maps_type &x_old, const Real &dt) const {

  // nlFunConstant b0("bmax", getValueFromMap(x_old, "b"));

  // Real Bint = _bdot.computeTimeIntegral(x, params, x_old, dt);

  // if (std::abs(Bint) < b0.computeValue(x, params, x_old))
  // return chain(std::copysign(1., Bint),
  //              _bdot.computeTimeIntegralParamGradient(x, params, x_old, dt));

  return _bdot.computeTimeIntegralParamGradient(x, params, x_old, dt);
  // else
  //   return b0.computeParamGradient(x, params, x_old);

  // io_maps_type b_values = {
  //     {"bmax", b0.computeValue(x, params, x_old)},
  //     {"bComp", _bdot.computeTimeIntegral(x, params, x_old, dt)}};
  //
  // io_maps_type b0_grad = b0.computeParamGradient(x, params, x_old);
  // io_maps_type b_grad =
  //     _bdot.computeTimeIntegralParamGradient(x, params, x_old, dt);
  //
  // map_of_io_maps_type b_values_grad = {{"bmax", &b0_grad}, {"bComp",
  // &b_grad}}; return _soft_min.computeParamGradient(b_values, b_values_grad,
  // params);
}
