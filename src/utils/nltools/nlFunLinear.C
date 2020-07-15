#pragma once

#include "nlFunLinear.h"

nlFunLinear::nlFunLinear(const std::string &fname, const std::string &var_name,
                         const Real m, const Real q)
    : nlFunBase(fname, requiredVar(), requiredParam()), _var_name(var_name),
      _m(m), _q(q) {}

Real nlFunLinear::computeValue(const io_maps_type &x,
                               const io_maps_type &params,
                               const io_maps_type &x_old) const {

  UNUSED(x_old);
  UNUSED(params);
  Real v = getValueFromMap(x, _var_name);
  return v * _m + _q;
}

nlFunLinear::io_maps_type
nlFunLinear::computeVarGradient(const io_maps_type &x,
                                const io_maps_type &params,
                                const io_maps_type &x_old) const {
  UNUSED(x);
  UNUSED(x_old);
  UNUSED(params);
  io_maps_type df_dvar = {{_var_name, _m}};

  return df_dvar;
}

nlFunLinear::io_maps_type
nlFunLinear::computeParamGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const {
  UNUSED(x);
  UNUSED(x_old);
  UNUSED(params);
  return _empty_map;
}
