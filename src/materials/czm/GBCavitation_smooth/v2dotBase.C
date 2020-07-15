#pragma once

#include "v2dotBase.h"

v2dotBase::v2dotBase(const std::string &fname, const Real &n, const Real &h,
                     const Real &theta_time_integration,
                     const equationTests &eqTests)
    : nlFunRateNoHistory(fname, theta_time_integration), _n(n), _h(h),
      _alpha(3. / (2. * n)), _eqTests(eqTests), _logit("logit", 1, 0.1, 0) {}

void v2dotBase::commonOperation(const io_maps_type &x,
                                const io_maps_type &params, Real &sVM, Real &sH,
                                Real &triax, Real &e, Real &a, Real &alpha,
                                Real &m, Real &beta) const {
  sVM = getValueFromMap(params, "sVM", "v2dotBase::params");
  sH = getValueFromMap(params, "sH", "v2dotBase::params");
  triax = 0;
  if (sVM != 0)
    triax = sH / sVM;
  e = getValueFromMap(params, "eqedotc", "v2dotBase::params");
  a = getValueFromMap(x, "a", "v2dotBase::x");
  m = std::copysign(1., sH);
  beta = betaFun(m);
  alpha = _alpha;
}

Real v2dotBase::gFun(const Real &m) const {
  Real g = 0;
  if (m < 0)
    g = 2. * _pi / (9. * std::pow(3., 0.5));
  else
    g = std::log(3.) - 2. / 3.;
  return g;
}
Real v2dotBase::betaFun(const Real &m) const {
  Real g = gFun(m);
  return (_n - 1) * (_n + g) / (_n * _n);
}
