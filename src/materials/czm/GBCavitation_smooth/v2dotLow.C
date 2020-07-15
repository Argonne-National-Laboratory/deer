#pragma once

#include "v2dotLow.h"

v2dotLow::v2dotLow(const std::string &fname, const Real &n, const Real &h,
                   const Real &theta_time_integration,
                   const equationTests &eqTests)
    : v2dotBase(fname, n, h, theta_time_integration, eqTests) {}

Real v2dotLow::computeValueLocal(const io_maps_type &x,
                                 const io_maps_type &params,
                                 const io_maps_type &x_old) const {
  UNUSED(x_old);
  Real sVM, sH, triax, e, a, alpha, m, beta;
  commonOperation(x, params, sVM, sH, triax, e, a, alpha, m, beta);

  Real f = 0;
  if (sVM != 0) {
    if (std::abs(triax) >= 1.) {
      f = 2. * e * std::pow(a, 3.) * _pi * _h * m *
          std::pow(alpha * std::abs(triax) + beta, _n);
    } else {
      f = 2. * e * std::pow(a, 3.) * _pi * _h * std::pow(alpha + beta, _n) *
          triax;
    }
  }
  if (std::isfinite(f))
    return f;
  else
    mooseError(_f_name + " value is not finite " + std::to_string(f));
  return 0;
}

v2dotLow::io_maps_type
v2dotLow::computeVarGradientLocal(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const {
  UNUSED(x_old);
  Real sVM, sH, triax, e, a, alpha, m, beta;
  commonOperation(x, params, sVM, sH, triax, e, a, alpha, m, beta);
  io_maps_type f_grad = _empty_map;

  if (sVM != 0) {
    if (std::abs(triax) >= 1.) {
      f_grad["a"] = 6. * e * std::pow(a, 2.) * _pi * _h * m *
                    std::pow(alpha * std::abs(triax) + beta, _n);

    } else {
      f_grad["a"] = 6. * e * std::pow(a, 2.) * _pi * _h *
                    std::pow(alpha + beta, _n) * triax;
    }
  }
  return f_grad;
}

v2dotLow::io_maps_type
v2dotLow::computeParamGradientLocal(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old) const {
  // UNUSED(x_old);
  // Real sVM, sH, triax, e, a, alpha, m, beta;
  // commonOperation(x, params, sVM, sH, triax, e, a, alpha, m, beta);
  // io_maps_type f_grad = _empty_map;
  // if (sVM != 0) {
  //   if (std::abs(triax) >= 1.) {
  //     f_grad["eqedotc"] = 2. * std::pow(a, 3.) * _pi * _h * m *
  //                         std::pow(alpha * std::abs(triax) + beta, _n);
  //
  //     f_grad["sH"] = 2 * _pi * std::pow(a, 3.) * alpha * e * _h * _n *
  //                    std::pow(alpha * std::abs(triax) + beta, _n) /
  //                    (alpha * std::abs(sH) + beta * std::abs(sVM));
  //
  //     f_grad["sVM"] =
  //         -2 * _pi * std::pow(a, 3.) * alpha * e * _h * m * _n *
  //         std::pow(alpha * std::abs(triax) + beta, _n) * std::abs(sH) /
  //         (std::abs(sVM) * (alpha * std::abs(sH) + beta * std::abs(sVM)));
  //
  //   } else {
  //
  //     f_grad["eqedotc"] =
  //         2. * std::pow(a, 3.) * _pi * _h * std::pow(alpha + beta, _n) *
  //         triax;
  //     f_grad["sH"] = 2. * e * std::pow(a, 3.) * _pi * _h *
  //                    std::pow(alpha + beta, _n) / sVM;
  //     f_grad["sVM"] = -2. * e * std::pow(a, 3.) * _pi * _h *
  //                     std::pow(alpha + beta, _n) * triax / sVM;
  //   }
  // }
  // return f_grad;
  return _empty_map;
}
