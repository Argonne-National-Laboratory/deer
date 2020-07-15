#pragma once

#include "f_aL.h"

f_aL::f_aL(const std::string &fname, const Real &D)
    : nlFunBase(fname, requiredVar(), requiredParam()), _D(D) {}

Real f_aL::computeValue(const io_maps_type &x, const io_maps_type &params,
                        const io_maps_type &x_old) const {
  UNUSED(x_old);

  Real sVM = getValueFromMap(params, "sVM");
  Real f = 0;
  if (sVM != 0) {
    Real a = getValueFromMap(x, "a");
    Real L = calc_L(params);
    f = std::pow(a / (a + 1.5 * L), 2.);
  }
  if (std::isfinite(f))
    return f;
  else
    mooseError(_f_name + "value is not finite " + std::to_string(f));
  return 0;
}

f_aL::io_maps_type f_aL::computeVarGradient(const io_maps_type &x,
                                            const io_maps_type &params,
                                            const io_maps_type &x_old) const {
  UNUSED(x_old);

  io_maps_type dfaL_dx = _empty_map;
  Real sVM = getValueFromMap(params, "sVM");
  if (sVM != 0) {
    Real a = getValueFromMap(x, "a");
    Real f = a * a;
    Real L = calc_L(params);
    Real g = std::pow(a + 1.5 * L, 2.);

    io_maps_type df_dx = {{"a", 2 * a}};
    io_maps_type dg_dx = {{"a", 2 * (a + 1.5 * L)}};

    dfaL_dx = f_divided_g_D(f, df_dx, g, dg_dx);
  }
  return dfaL_dx;
}

f_aL::io_maps_type f_aL::computeParamGradient(const io_maps_type &x,
                                              const io_maps_type &params,
                                              const io_maps_type &x_old) const {
  UNUSED(x_old);
  io_maps_type df_dx = _empty_map;
  return df_dx;
}

Real f_aL::calc_L(const io_maps_type &params) const {
  Real e = getValueFromMap(params, "eqedotc");
  Real sVM = getValueFromMap(params, "sVM");

  return std::pow(_D * sVM / e, 1. / 3.);
}
