#pragma once

#include "logistic.h"

logistic::logistic(const std::string &fname, const Real &L, const Real &k,
                   const Real &x0)
    : nlFunBase(fname), _L(L), _k(k), _x0(x0) {}

Real logistic::computeLogisticValue(const Real &x) const {

  return _L / (1 + std::exp(-_k * (x - _x0)));
}

logistic::io_maps_type logistic::computeLogisticGradient(const Real &x) const {

  io_maps_type dlogistic_dx;
  Real e = std::exp(-_k * (x - _x0));
  dlogistic_dx["x"] = _L * _k * e / std::pow(1. + e, 2.);
  return dlogistic_dx;
}
