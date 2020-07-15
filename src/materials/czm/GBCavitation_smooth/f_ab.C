#pragma once

#include "f_ab.h"

f_ab::f_ab(const std::string &fname)
    : nlFunBase(fname, requiredVar(), requiredParam()) {}

Real f_ab::computeValue(const io_maps_type &x, const io_maps_type &params,
                        const io_maps_type &x_old) const {
  UNUSED(x_old);
  UNUSED(params);
  Real f = std::pow(getValueFromMap(x, "a") / getValueFromMap(x, "b"), 2.);

  if (f == 1)
    return f * (1. - 1e-6);

  if (std::isfinite(f))
    return f;

  else
    mooseError(_f_name + "value is not finite " + std::to_string(f));
  return 0;
}

f_ab::io_maps_type f_ab::computeVarGradient(const io_maps_type &x,
                                            const io_maps_type &params,
                                            const io_maps_type &x_old) const {
  UNUSED(x_old);
  UNUSED(params);

  Real a = getValueFromMap(x, "a");
  io_maps_type df_dx = {{"a", 2. * a}};
  Real b = getValueFromMap(x, "b");
  io_maps_type dg_dx = {{"b", 2. * b}};
  return f_divided_g_D(a * a, df_dx, b * b, dg_dx);
}
