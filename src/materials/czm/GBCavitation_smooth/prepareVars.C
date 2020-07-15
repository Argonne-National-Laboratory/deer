#pragma once

#include "prepareVars.h"

prepareVars::prepareVars(const std::string &fname, const Real &a0)
    : nlFunBase(fname, requiredVar(), requiredParam()), _a0(a0) {}

prepareVars::io_maps_type
prepareVars::computeNewX(const io_maps_type &x_NL,
                         const io_maps_type &x_old) const {
  io_maps_type newx = x_NL;
  // Real aold = getValueFromMap(x_old, "a");
  // Real bold = getValueFromMap(x_old, "b");
  //
  // newx["a"] = std::abs(getValueFromMap(x_NL, "a")) * aold + _a0;
  // newx["b"] = -1. * std::abs(getValueFromMap(x_NL, "b")) * bold + bold;
  newx["a"] = std::abs(getValueFromMap(x_NL, "a"));
  newx["b"] = std::abs(getValueFromMap(x_NL, "b"));
  return newx;
}

prepareVars::io_maps_type
prepareVars::computeNewXVarGradient(const io_maps_type &xNL,
                                    const io_maps_type &x_old) const {

  Real aold = getValueFromMap(x_old, "a");
  Real bold = getValueFromMap(x_old, "b");

  io_maps_type dx_dxNL = xNL;
  for (auto &it : dx_dxNL)
    it.second = 1.;
  // dx_dxNL["a"] = std::copysign(1., getValueFromMap(xNL, "a")) * aold;
  // dx_dxNL["b"] = -1. * std::copysign(1., getValueFromMap(xNL, "b")) * bold;
  dx_dxNL["a"] = std::copysign(1., getValueFromMap(xNL, "a"));
  dx_dxNL["b"] = std::copysign(1., getValueFromMap(xNL, "b"));
  return dx_dxNL;
}

prepareVars::io_maps_type
prepareVars::XNLFromComputedX(const io_maps_type &comp_x,
                              const io_maps_type &x_old) const {
  io_maps_type newNLx = comp_x;
  // Real aold = getValueFromMap(x_old, "a");
  // Real bold = getValueFromMap(x_old, "b");
  //
  // newNLx["a"] = (getValueFromMap(comp_x, "a") - _a0) / aold;
  // newNLx["b"] = (getValueFromMap(comp_x, "b") - bold) / bold;
  return newNLx;
}
