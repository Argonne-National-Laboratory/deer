#pragma once

#include "pNormSoftMax.h"

pNormSoftMax::pNormSoftMax(const std::string &fname, const Real &p)
    : nlFunBase(fname), _p(p) {}

Real pNormSoftMax::computeSoftMaxValue(const io_maps_type &x) const {

  Real max = 0;
  for (auto &it : x)
    max += std::pow(it.second, _p);
  return std::pow(max, 1. / _p);
}

pNormSoftMax::io_maps_type
pNormSoftMax::computeSoftMaxGradient(const io_maps_type &x) const {

  io_maps_type max_grad;
  Real sum_power = 0;
  for (auto &it : x)
    sum_power += std::pow(it.second, _p);

  for (auto &it : x)
    max_grad[it.first] =
        std::pow(sum_power, 1. / _p - 1.) * std::pow(it.second, _p - 1.);

  return max_grad;
}
