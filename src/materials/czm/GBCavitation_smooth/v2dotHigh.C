#pragma once

#include "v2dotHigh.h"

v2dotHigh::v2dotHigh(const std::string &fname, const Real &n, const Real &h,
                     const Real &theta_time_integration,
                     const equationTests &eqTests)
    : v2dotBase(fname, n, h, theta_time_integration, eqTests) {}

void v2dotHigh::commonOperation(const io_maps_type &x,
                                const io_maps_type &params, Real &sVM, Real &sH,
                                Real &triax, Real &e, Real &a, Real &alpha,
                                Real &m, Real &beta, Real &b, Real &Kab_pow,
                                Real &Kab_pow_1, Real &num_arg_n,
                                Real &pow_arg_n) const {

  // Real sVM, sH, triax, e, a, alpha, m, beta;
  v2dotBase::commonOperation(x, params, sVM, sH, triax, e, a, alpha, m, beta);

  // Real b, Kab_pow, Kab_pow_1, num_arg_n, pow_arg_n;
  b = getValueFromMap(x, "b", "v2dotBase::x");
  Kab_pow = std::pow(0.87 * a / b, 3. / _n);
  Kab_pow_1 = Kab_pow - 1;
  if (std::abs(triax) > 1)
    num_arg_n = alpha * _n * std::abs(sH) + m * std::abs(sVM);
  else
    num_arg_n = alpha * _n + m;

  pow_arg_n = 0;
  if (sVM != 0)
    pow_arg_n = std::pow((-num_arg_n) / (_n * Kab_pow_1 * std::abs(sVM)), _n);
}

Real v2dotHigh::computeValueLocal(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const {
  UNUSED(x_old);
  Real sVM, sH, triax, e, a, alpha, m, beta;
  v2dotBase::commonOperation(x, params, sVM, sH, triax, e, a, alpha, m, beta);

  Real num;
  Real b = getValueFromMap(x, "b");
  Real f = 0;
  Real den = 1. - std::pow(0.87 * a / b, 3. / _n);
  if (sVM != 0) {
    if (std::abs(triax) >= 1.) {
      num = alpha * std::abs(triax) + m / _n;
      f = 2. * e * std::pow(a, 3.) * _pi * _h * m * std::pow(num / den, _n);
    } else {
      num = alpha + m / _n;
      f = 2. * e * std::pow(a, 3.) * _pi * _h * std::pow(num / den, _n) * triax;
    }
  }
  if (std::isfinite(f))
    return f;
  else {
    mooseError(_f_name + " value is not finite " + std::to_string(f));
    return 0;
  }
}

v2dotHigh::io_maps_type
v2dotHigh::computeVarGradientLocal(const io_maps_type &x,
                                   const io_maps_type &params,
                                   const io_maps_type &x_old) const {
  UNUSED(x_old);

  io_maps_type f_grad = _empty_map;
  Real sVM, sH, triax, e, a, alpha, m, beta;
  v2dotBase::commonOperation(x, params, sVM, sH, triax, e, a, alpha, m, beta);
  if (sVM != 0) {

    Real b = getValueFromMap(x, "b");
    Real P, num, dP_da, dP_db;
    Real den = 1. - std::pow(0.87 * a / b, 3. / _n);
    io_maps_type dden_dx = _empty_map;
    io_maps_type dnum_dx = _empty_map;
    Real prefactor = 2. * e * std::pow(a, 3.) * _pi * _h * m;
    io_maps_type dprefactor_dx = {
        {"a", 6. * e * std::pow(a, 2.) * _pi * _h * m}};

    dden_dx["a"] = -3. / _n * std::pow(0.87 * a / b, 3. / _n - 1.) * 0.87 / b;
    dden_dx["b"] =
        3. / _n * std::pow(0.87 * a / b, 3. / _n - 1.) * 0.87 * a / (b * b);
    if (std::abs(triax) >= 1.)
      num = alpha * std::abs(triax) + m / _n;
    else
      num = alpha + m / _n;

    P = std::pow(num / den, _n);
    /* numerator is constant, no derivatives we can compute the deriuvatives of
     * C^n*f^-n */
    io_maps_type dP_dx =
        chain(std::pow(num, _n), f_power_n_D(den, dden_dx, -_n));

    f_grad = f_times_g_D(prefactor, dprefactor_dx, P, dP_dx);
  }
  return f_grad;
}

v2dotHigh::io_maps_type
v2dotHigh::computeParamGradientLocal(const io_maps_type &x,
                                     const io_maps_type &params,
                                     const io_maps_type &x_old) const {
  // UNUSED(x_old);
  // Real sVM, sH, triax, e, a, alpha, m, beta, b, Kab_pow, Kab_pow_1,
  // num_arg_n,
  //     pow_arg_n;
  // commonOperation(x, params, sVM, sH, triax, e, a, alpha, m, beta, b,
  // Kab_pow,
  //                 Kab_pow_1, num_arg_n, pow_arg_n);
  // io_maps_type f_grad = _empty_map;
  // if (sVM != 0) {
  //   if (std::abs(triax) >= 1.) {
  //     f_grad["eqedotc"] = 2. * _pi * std::pow(a, 3.) * _h * m * pow_arg_n;
  //     f_grad["sH"] = 2. * _pi * std::pow(a, 3.) * alpha * e * _h *
  //                    std::pow(_n, 2.) * pow_arg_n / num_arg_n; /* times m^2*/
  //     f_grad["sVM"] = -2. * _pi * std::pow(a, 3.) * alpha * e * _h * m *
  //                     std::pow(_n, 2.) * pow_arg_n * std::abs(sH) /
  //                     (sVM * num_arg_n);
  //   } else {
  //     f_grad["eqedotc"] = 2. * _pi * std::pow(a, 3.) * _h * triax *
  //     pow_arg_n; f_grad["sH"] = 2. * _pi * std::pow(a, 3.) * e * _h *
  //     pow_arg_n / sVM; f_grad["sVM"] =
  //         -2. * _pi * std::pow(a, 3.) * e * _h * triax * pow_arg_n / sVM;
  //   }
  // }
  // return f_grad;
  return _empty_map;
}
