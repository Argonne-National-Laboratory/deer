#pragma once

#include "bdot.h"

bdot::bdot(const std::string &f_name, const Real &beta, const Real &b0,
           const Real &bsat, const Real &Sthr, const Real &sigma0,
           const Real &FN, const equationTests &eqTest,
           const Real &theta_time_integration)
    : nlFunRate(f_name, "b", theta_time_integration), _beta(beta), _b0(b0),
      _bsat(bsat), _Sthr(Sthr), _sigma0(sigma0), _FN(FN), _eqTest(eqTest),
      _soft_max_min("soft_max_min", 100), _soft_max("soft_max", 35) {}

Real bdot::computeValue(const io_maps_type &x, const io_maps_type &params,
                        const io_maps_type &x_old) const {

  if (!_eqTest.nucleationIsActive(x, params, x_old))
    return 0;

  Real TN = getValueFromMap(x, "Tn", "bdot::computeValue.x");
  Real sTN = _soft_max.computeSoftMaxValue(_TNmin, TN);

  Real b = getValueFromMap(x, "b", "bdot::computeValue.x");
  Real eqedotc =
      getValueFromMap(params, "eqedotc", "bdot::computeValue.params");
  Real b_dot =
      -_pi * std::pow(b, 3.) * _FN * std::pow(sTN / _sigma0, _beta) * eqedotc;

  return b_dot;
}

bdot::io_maps_type bdot::computeVarGradient(const io_maps_type &x,
                                            const io_maps_type &params,
                                            const io_maps_type &x_old) const {
  io_maps_type f_grad = _empty_map;
  if (!_eqTest.nucleationIsActive(x, params, x_old))
    return _empty_map;

  Real TN = getValueFromMap(x, "Tn", "bdot::computeValue.x");
  Real sTN = _soft_max.computeSoftMaxValue(_TNmin, TN);
  Real dsTN_dTN =
      _soft_max.computeSoftMaxGradient(_TNmin, TN).find("b")->second;

  Real b = getValueFromMap(x, "b", "bdot::computeVarGradient.x");
  Real eqedotc =
      getValueFromMap(params, "eqedotc", "bdot::computeVarGradient.params");
  f_grad["b"] = -3. * _pi * std::pow(b, 2.) * _FN *
                std::pow(sTN / _sigma0, _beta) * eqedotc;
  f_grad["Tn"] = -_beta * _pi * std::pow(b, 3.) * _FN *
                 std::pow(sTN / _sigma0, _beta - 1.) * eqedotc / _sigma0 *
                 dsTN_dTN;

  return f_grad;
}

bdot::io_maps_type bdot::computeParamGradient(const io_maps_type &x,
                                              const io_maps_type &params,
                                              const io_maps_type &x_old) const {

  return _empty_map;
}
