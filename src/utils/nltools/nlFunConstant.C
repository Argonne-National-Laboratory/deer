#pragma once

#include "nlFunConstant.h"

nlFunConstant::nlFunConstant(const std::string &fname, const Real &C)
    : nlFunBase(fname, requiredVar(), requiredParam()), _C(C) {
  _f_value = _C;
  _f_var_gradient = _empty_map;
  _f_param_gradient = _empty_map;
}

Real nlFunConstant::computeValue(const io_maps_type &x,
                                 const io_maps_type &params,
                                 const io_maps_type &x_old) const {
  UNUSED(x);
  UNUSED(x_old);
  UNUSED(params);
  return _C;
}

nlFunConstant::io_maps_type
nlFunConstant::computeVarGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const {
  UNUSED(x);
  UNUSED(x_old);
  UNUSED(params);
  return _empty_map;
}

nlFunConstant::io_maps_type
nlFunConstant::computeParamGradient(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old) const {
  UNUSED(x);
  UNUSED(x_old);
  UNUSED(params);
  return _empty_map;
}
