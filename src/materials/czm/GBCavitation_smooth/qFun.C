#pragma once

#include "f_Low.h"
#include "f_ab.h"
#include "nlFunConstant.h"
#include "qFun.h"

qFun::qFun(const std::string &fname, const nlFunBase &f)
    : nlFunBase(fname, requiredVar(), requiredParam()), _f(f) {}

Real qFun::computeValue(const io_maps_type &x, const io_maps_type &params,
                        const io_maps_type &x_old) const {
  UNUSED(x_old);
  Real f = _f.computeValue(x, params, x_old);
  Real q = 2 * std::log(1. / f) - (1. - f) * (3. - f);
  if (q == 0)
    mooseError(_f_name + " q==0!! f = " + std::to_string(f));
  return q;
}

qFun::io_maps_type qFun::computeVarGradient(const io_maps_type &x,
                                            const io_maps_type &params,
                                            const io_maps_type &x_old) const {
  UNUSED(x_old);
  Real dq_df = dqdf(x, params, x_old);
  io_maps_type dq_dx = chain(dq_df, _f.computeVarGradient(x, params, x_old));
  return dq_dx;
}

qFun::io_maps_type qFun::computeParamGradient(const io_maps_type &x,
                                              const io_maps_type &params,
                                              const io_maps_type &x_old) const {
  UNUSED(x_old);
  return _empty_map;
}

Real qFun::dqdf(const io_maps_type &x, const io_maps_type &params,
                const io_maps_type &x_old) const {
  Real f = _f.computeValue(x, params, x_old);
  return -2. * f + 4. - 2. / f;
}
