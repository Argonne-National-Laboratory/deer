#pragma once

#include "aComp.h"

aComp::aComp(const std::string &fname, const adot &aFun, const Real &amin,
             const Real &theta_time_integration)
    : nlFunRate(fname, "a", theta_time_integration), _adot(aFun),
      _a0Fun("amin", amin), _soft_max_min("soft_maxmin", 100),
      _soft_max("soft_maxmin", 35) {}

Real aComp::computeTimeIntegral(const io_maps_type &x,
                                const io_maps_type &params,
                                const io_maps_type &x_old,
                                const Real &dt) const {

  // io_maps_type a_values = {
  //     {"amin", _a0Fun.computeValue(x, params, x_old)},
  //     {"aComp", _adot.computeTimeIntegral(x, params, x_old, dt)}};

  Real f = _soft_max.computeSoftMaxValue(
      _a0Fun.computeValue(x, params, x_old),
      _adot.computeTimeIntegral(x, params, x_old, dt));
  // Real f = _soft_max.computeValue(a_values, params, x_old);

  if (std::isfinite(f))
    return f;
  else
    mooseError(_f_name + "value is not finite " + std::to_string(f));
  return 0;
}

aComp::io_maps_type aComp::computeTimeIntegralVarGradient(
    const io_maps_type &x, const io_maps_type &params,
    const io_maps_type &x_old, const Real &dt) const {

  io_maps_type sm_grad = _soft_max.computeSoftMaxGradient(
      _a0Fun.computeValue(x, params, x_old),
      _adot.computeTimeIntegral(x, params, x_old, dt));
  // Real amin = _a0Fun.computeValue(x, params, x_old);
  // Real a = _adot.computeTimeIntegral(x, params, x_old, dt);

  // io_maps_type a_values = {
  //     {"amin", _a0Fun.computeValue(x, params, x_old)},
  //     {"aComp", _adot.computeTimeIntegral(x, params, x_old, dt)}};
  //
  // io_maps_type a0_grad = _a0Fun.computeVarGradient(x, params, x_old);
  // io_maps_type a_grad =
  //     _adot.computeTimeIntegralVarGradient(x, params, x_old, dt);
  // map_of_io_maps_type a_values_grad = {{"amin", &a0_grad}, {"aComp",
  // &a_grad}};

  // return _soft_max.computeVarGradient(a_values, a_values_grad, params);
  return chain(sm_grad.find("b")->second,
               _adot.computeTimeIntegralVarGradient(x, params, x_old, dt));
}

aComp::io_maps_type aComp::computeTimeIntegralParamGradient(
    const io_maps_type &x, const io_maps_type &params,
    const io_maps_type &x_old, const Real &dt) const {
  // io_maps_type a_values = {
  //     {"amin", _a0Fun.computeValue(x, params, x_old)},
  //     {"aComp", _adot.computeTimeIntegral(x, params, x_old, dt)}};
  //
  // io_maps_type a0_grad = _a0Fun.computeParamGradient(x, params, x_old);
  // io_maps_type a_grad =
  //     _adot.computeTimeIntegralParamGradient(x, params, x_old, dt);
  //
  // map_of_io_maps_type a_values_grad = {{"amin", &a0_grad}, {"aComp",
  // &a_grad}}; return _soft_max.computeParamGradient(a_values, a_values_grad,
  // params);

  return _adot.computeTimeIntegralParamGradient(x, params, x_old, dt);
}
