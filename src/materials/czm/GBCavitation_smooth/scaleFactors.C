#pragma once

#include "scaleFactors.h"

scaleFactors::scaleFactors(const std::string &fname)
    : nlFunBase(fname, requiredVar(), requiredParam()) {}

scaleFactors::io_maps_type
scaleFactors::computeVarScaleFactor(const io_maps_type &x_old) const {
  io_maps_type scale_factor = _empty_map;

  scale_factor["a"] = getValueFromMap(x_old, "a");
  scale_factor["b"] = getValueFromMap(x_old, "b");
  scale_factor["Tn"] = std::abs(getValueFromMap(x_old, "Tn"));
  // scale_factor["a"] = 0.001;
  // scale_factor["b"] = 0.001;
  // scale_factor["Tn"] = 10;
  scale_factor["Ts1"] = 1;
  scale_factor["Ts2"] = 1;

  for (auto &it : scale_factor) {
    if (it.second == 0)
      it.second = 1;
    if (it.first == "Tn")
      if (std::abs(it.second) < 10.)
        it.second = 1.;
  }
  return scale_factor;
}
