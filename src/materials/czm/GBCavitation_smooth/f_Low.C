#pragma once

#include "f_Low.h"

f_Low::f_Low(const std::string &fname, const f_aL &faL, const f_ab &fab)
    : nlFunBase(fname, requiredVar(), requiredParam()), _faL(faL), _fab(fab),
      _soft_max_min("soft_max_min", 100,
                    /*max_of_abs =*/false,
                    /*const bool &max_of_abs_signed = */ false,
                    /*const bool &compute_min = */ false),
      _soft_max("soft_max") {}

Real f_Low::computeValue(const io_maps_type &x, const io_maps_type &params,
                         const io_maps_type &x_old) const {
  UNUSED(x_old);

  Real fal = _faL.computeValue(x, params, x_old);
  Real fab = _fab.computeValue(x, params, x_old);

  return _soft_max.computeSoftMaxValue(fal, fab);
}

f_Low::io_maps_type f_Low::computeVarGradient(const io_maps_type &x,
                                              const io_maps_type &params,
                                              const io_maps_type &x_old) const {

  Real fal = _faL.computeValue(x, params, x_old);
  Real fab = _fab.computeValue(x, params, x_old);
  io_maps_type dfaL_dx = _faL.computeVarGradient(x, params, x_old);
  io_maps_type dfab_dx = _fab.computeVarGradient(x, params, x_old);

  io_maps_type dmax_dfi = _soft_max.computeSoftMaxGradient(fal, fab);

  Real dmax_dfal = dmax_dfi.find("a")->second;
  Real dmax_dfab = dmax_dfi.find("b")->second;
  dfaL_dx = chain(dmax_dfal, dfaL_dx);
  dfab_dx = chain(dmax_dfab, dfab_dx);
  return sumD(dfaL_dx, dfab_dx);
}

f_Low::io_maps_type
f_Low::computeParamGradient(const io_maps_type &x, const io_maps_type &params,
                            const io_maps_type &x_old) const {
  UNUSED(x_old);

  return _empty_map;
}
