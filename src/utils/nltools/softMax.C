#pragma once

#include "softMax.h"

softMax::softMax(const std::string &fname, const Real &f)
    : nlFunBase(fname), _f(f) {}

Real softMax::computeSoftMaxValue(const Real &a, const Real &b) const {

  Real s = std::abs(a) + std::abs(b);
  return std::log(std::exp(a / s * _f) + std::exp(b / s * _f)) / _f * s;
}

softMax::io_maps_type softMax::computeSoftMaxGradient(const Real &a,
                                                      const Real &b) const {

  io_maps_type dsmax_dx;
  Real a_abs = std::abs(a);
  Real a_sign = std::copysign(1., a);
  Real b_abs = std::abs(b);
  Real b_sign = std::copysign(1., b);
  Real s = std::abs(a) + std::abs(b);

  Real f_s = _f / s;
  Real f_s2 = _f / (s * s);
  Real exp_a = std::exp(a * f_s);
  Real exp_b = std::exp(b * f_s);
  Real f_sumexp = _f * (exp_a + exp_b);

  Real term1 = std::log(exp_a + exp_b) / _f;

  /*dmax_da*/
  Real term1_a = term1 * a_sign;
  Real term2_a =
      (exp_a * (f_s - a * a_sign * f_s2) - (b * exp_b * a_sign * f_s2)) *
      (a_abs + b_abs) / f_sumexp;
  dsmax_dx["a"] = term1_a + term2_a;

  /*dmax_db*/
  Real term1_b = term1 * b_sign;
  Real term2_b =
      (exp_b * (f_s - b * b_sign * f_s2) - (a * exp_a * b_sign * f_s2)) *
      (a_abs + b_abs) / f_sumexp;
  dsmax_dx["b"] = term1_b + term2_b;

  return dsmax_dx;
}
