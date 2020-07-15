#pragma once

#include "softMaxMin.h"

softMaxMin::softMaxMin(const std::string &fname, const Real &f /*= 50*/,
                       const bool &max_of_abs /*= false*/,
                       const bool &max_of_abs_signed /*= false*/,
                       const bool &compute_min /*= false*/)
    : nlFunBase(fname), _f(f), _max_of_abs(max_of_abs),
      _max_of_abs_signed(max_of_abs_signed), _compute_min(compute_min) {

  // check inputs
  if (_max_of_abs && _max_of_abs_signed)
    mooseError("softMaxMin: only one between max_of_abs and "
               "max_of_abs_signed can be true");
}

Real softMaxMin::computeValue(const io_maps_type &x, const io_maps_type &params,
                              const io_maps_type &x_old) const {

  UNUSED(x_old);
  UNUSED(params);
  io_maps_type _x2use, _sign_x2use;
  Real _sign_xref, _c, _k, _arg_log;
  updateX2useAndParam(x, _x2use, _sign_x2use, _sign_xref);
  updateParam(_x2use, _c, _k, _arg_log);
  Real smax = 0;
  smax = std::log(_arg_log) / _k + _c;
  if (_max_of_abs_signed)
    smax *= _sign_xref;
  if (_compute_min)
    smax *= -1;

  return smax;
}

softMaxMin::io_maps_type
softMaxMin::computeVarGradient(const io_maps_type &x,
                               const io_maps_type &params,
                               const io_maps_type &x_old) const {
  UNUSED(x_old);
  UNUSED(params);
  io_maps_type _x2use, _sign_x2use;
  Real _sign_xref, _c, _k, _arg_log;
  updateX2useAndParam(x, _x2use, _sign_x2use, _sign_xref);
  updateParam(_x2use, _c, _k, _arg_log);

  io_maps_type dsmax_dx = x;
  io_maps_type grad_max;
  grad_max.clear();
  std::unordered_map<std::string, io_maps_type> dfi_dxi_temp;

  Real sum_ekx = 0;
  Real sum_czekxc = 0;
  for (const auto &it : _x2use) {
    Real xv = getValueFromMap(_x2use, it.first, "_x2use");
    sum_ekx += std::exp(xv * _k);
    sum_czekxc += (xv - _c) * std::exp(_k * (xv - _c));
  }
  Real dM_dk = sum_czekxc / (_k * _arg_log) - std::log(_arg_log) / (_k * _k);

  for (const auto &it : _x2use) {
    Real xv = getValueFromMap(_x2use, it.first, "_x2use");
    Real xs = getValueFromMap(_sign_x2use, it.first, "_sign_x2use");
    Real dk_dxi = -_f * std::copysign(1., xv) / (_c * _c);
    Real dM_dxi = std::exp(xv * _k) / (sum_ekx);
    Real dM_dxi_total = dM_dxi + dM_dk * dk_dxi;
    if (_max_of_abs || _max_of_abs_signed)
      dM_dxi_total *= xs;
    if (_max_of_abs_signed)
      dM_dxi_total *= _sign_xref;

    setValueInMap(dsmax_dx, it.first, dM_dxi_total, "dsmax_dx");
  }

  return dsmax_dx;
}

softMaxMin::io_maps_type
softMaxMin::computeVarGradient(const io_maps_type &x,
                               const map_of_io_maps_type &x_grad,
                               const io_maps_type &params) const {
  io_maps_type total_grad;

  io_maps_type max_grad = computeVarGradient(x, params, _empty_map);

  for (const auto &it : x_grad) {
    Real dMax_dfi = getValueFromMap(max_grad, it.first,
                                    "softMaxMin::computeGradient::max_grad");
    io_maps_type dfi_dxi = chain(dMax_dfi, *it.second);
    total_grad = sumD(total_grad, dfi_dxi);
  }
  return total_grad;
}

softMaxMin::io_maps_type
softMaxMin::computeParamGradient(const io_maps_type &x,
                                 const map_of_io_maps_type &x_grad,
                                 const io_maps_type &params) const {
  io_maps_type total_grad;

  io_maps_type max_grad = computeVarGradient(x, params, _empty_map);

  for (const auto &it : x_grad) {
    Real dMax_dfi = getValueFromMap(max_grad, it.first,
                                    "softMaxMin::computeGradient::max_grad");
    io_maps_type dfi_dxi = chain(dMax_dfi, *it.second);
    total_grad = sumD(total_grad, dfi_dxi);
  }
  return total_grad;
}

void softMaxMin::updateParam(const io_maps_type &_x2use, Real &_c, Real &_k,
                             Real &_arg_log) const {
  // compute c
  _c = 0;
  for (const auto &it : _x2use)
    _c += std::abs(getValueFromMap(_x2use, it.first));
  _k = 1 / _c * _f;
  _arg_log = 0;
  // compute argument of the log
  for (const auto &it : _x2use)
    _arg_log += std::exp((getValueFromMap(_x2use, it.first) - _c) * _k);
}

void softMaxMin::updateX2useAndParam(const io_maps_type &x,
                                     io_maps_type &_x2use,
                                     io_maps_type &_sign_x2use,
                                     Real &_sign_xref) const {
  _x2use = x;
  _sign_x2use = x;
  for (auto &it : x)
    setValueInMap(_sign_x2use, it.first, std::copysign(1, it.second),
                  "_sign_x2use");

  _sign_xref = 1;
  Real xref = std::abs(x.begin()->second);

  if (_max_of_abs_signed) {
    for (const auto &it : x) {
      Real xv = getValueFromMap(x, it.first);
      Real xv_abs = std::abs(xv);
      Real xs = std::copysign(1, xv);
      if ((_compute_min && (xv_abs <= xref)) ||
          (!_compute_min && (xv_abs >= xref))) {
        xref = xv_abs;
        _sign_xref = xs;
      }
    }
  }

  if (_max_of_abs || _max_of_abs_signed) {
    for (const auto &it : x) {
      Real xv = getValueFromMap(x, it.first);
      setValueInMap(_x2use, it.first, std::abs(xv), "_x2use");
      // setValueInMap(_sign_x2use, it.first, std::copysign(1, xv),
      // "_sign_x2use");
    }
  }

  if (_compute_min) {
    for (const auto &it : _x2use) {
      Real xv = getValueFromMap(_x2use, it.first);
      setValueInMap(_x2use, it.first, -1 * xv, "_x2use");
    }
  }
}

// static void softMaxMin(const io_maps_type &x,
//                        Real &smax,
//                        io_maps_type &dsmax_dx,
//                        const bool &compute_gradient = true, const Real &f =
//                        50, const bool &max_of_abs = false, const bool
//                        &max_of_abs_signed = false, const bool &compute_min =
//                        false) {
//
//   if (max_of_abs && max_of_abs_signed)
//     mooseError("smoothMaxExp: only one between max_of_abs and "
//                "max_of_abs_signed can be true");
//
//   smax = 0;
//   dsmax_dx.clear();
//   dsmax_dx.insert(x.begin(), x.end());
//
//   io_maps_type x2use = x;
//   io_maps_type sign_x2use = x;
//
//   Real xref = std::abs(x.begin()->second);
//   Real sign_xref = 0;
//
//   if (max_of_abs_signed) {
//     for (const auto &it : x) {
//       Real xv = getValueFromMap(x, it.first);
//       Real xv_abs = std::abs(xv);
//       Real xs = std::copysign(1, xv);
//       if ((compute_min && (xv_abs <= xref)) ||
//           (!compute_min && (xv_abs >= xref))) {
//         xref = xv_abs;
//         sign_xref = xv;
//       }
//     }
//   }
//
//   if (max_of_abs || max_of_abs_signed) {
//     for (const auto &it : x) {
//       Real xv = getValueFromMap(x, it.first);
//       setValueInMap(x2use, it.first, std::abs(xv), "x2use");
//       setValueInMap(sign_x2use, it.first, std::copysign(1, xv),
//       "sign_x2use");
//     }
//   }
//
//   if (compute_min) {
//     for (const auto &it : x2use) {
//       Real xv = getValueFromMap(x2use, it.first);
//       setValueInMap(x2use, it.first, -1 * xv, "x2use");
//     }
//   }
//
//   Real arg_log = 0;
//   Real c = 0;
//   Real k = 0;
//
//   // compute c
//   for (const auto &it : x2use)
//     c += std::abs(getValueFromMap(x2use, it.first));
//   k = 1 / c * f;
//   // compute argument of the log
//   for (const auto &it : x2use)
//     arg_log += exp((getValueFromMap(x2use, it.first) - c) * k);
//
//   smax = std::log(arg_log) / k + c;
//   if (max_of_abs_signed == 1)
//     smax *= sign_xref;
//   if (compute_min == 1)
//     smax *= -1;
//
//   if (compute_gradient) {
//     Real sum_ekx = 0;
//     Real sum_czekxc = 0;
//     for (const auto &it : x2use) {
//       Real xv = getValueFromMap(x2use, it.first);
//       sum_ekx += std::exp(xv * k);
//       sum_czekxc += (xv - c) * std::exp(k * (xv - c));
//     }
//     Real dM_dk = sum_czekxc / (k * arg_log) - std::log(arg_log) / (k * k);
//
//     for (const auto &it : x2use) {
//       Real xv = getValueFromMap(x2use, it.first);
//       Real xs = getValueFromMap(sign_x2use, it.first);
//       Real dk_dxi = -f * xs / (c * c);
//       Real dM_dxi = std::exp(xv * k) / (sum_ekx);
//       Real dM_dxi_total = dM_dxi + dM_dk * dk_dxi;
//       if (max_of_abs || max_of_abs_signed) {
//         dM_dxi_total *= xs;
//         if (max_of_abs_signed)
//           dM_dxi_total *= sign_xref;
//       }
//       setValueInMap(dsmax_dx, it.first, dM_dxi_total, "dsmax_dx");
//     }
//   }
// }

// static void softMaxMinNumericalGradient(
//     const io_maps_type &x,
//     io_maps_type &dmax_dx_numerical,
//     const Real &f = 50, const bool &max_of_abs = false,
//     const bool &max_of_abs_signed = false, const bool &compute_min = false) {
//
//   if (max_of_abs && max_of_abs_signed)
//     mooseError("smoothMaxExp: only one between max_of_abs and "
//                "max_of_abs_signed can be true");
//
//   Real soft_max0 = 0;
//
//   io_maps_type dmax_df_temp;
//
//   softMaxMin(x, soft_max0, dmax_dx_numerical,
//              /*compute_gradient = */ false, f, max_of_abs, max_of_abs_signed,
//              compute_min);
//
//   dmax_dx_numerical.clear();
//   dmax_dx_numerical.insert(x.begin(), x.end());
//   io_maps_type x_temp;
//
//   Real x_plus_dx, dx, soft_max_temp, x_val_temp;
//   Real dM_df;
//   for (const auto &it : x) {
//     soft_max_temp = 0;
//     x_temp.clear();
//     x_temp = x;
//     x_val_temp = getValueFromMap(x_temp, it.first, "x_temp");
//     optimalDx(x_val_temp, x_plus_dx, dx);
//     setValueInMap(x_temp, it.first, x_plus_dx, "x_temp");
//     softMaxMin(x_temp, soft_max_temp, dmax_df_temp,
//                /*compute_gradient = */ false, f, max_of_abs,
//                max_of_abs_signed, compute_min);
//     dM_df = numericDiff(soft_max0, soft_max_temp, dx);
//     setValueInMap(dmax_dx_numerical, it.first, dM_df, "dmax_dx_numerical");
//   }
// }
//
// static void
// checksoftMaxMinGradient(const io_maps_type &x,
//                         const Real &f = 50, const bool &max_of_abs = false,
//                         const bool &max_of_abs_signed = false,
//                         const bool &compute_min = false) {
//
//   if (max_of_abs && max_of_abs_signed)
//     mooseError("checksoftMaxMinGradient: only one between max_of_abs and "
//                "max_of_abs_signed can be true");
//
//   Real smax_0;
//   io_maps_type dmax_df, dmax_df_numerical;
//   softMaxMin(x, smax_0, dmax_df,
//              /*compute_gradient = */ true, f, max_of_abs, max_of_abs_signed,
//              compute_min);
//
//   softMaxMinNumericalGradient(x, dmax_df_numerical, f, max_of_abs,
//                               max_of_abs_signed, compute_min);
//
//   for (const auto &it : dmax_df) {
//     Real num_der =
//         getValueFromMap(dmax_df_numerical, it.first, "dmax_df_numerical");
//     Real der = nltoolsns::getValueFromMap(dmax_df, it.first, "dmax_df");
//     if (std::copysign(1, num_der) != std::copysign(1, der))
//       Moose::out << "checksoftMaxMinGradient the sign of the numerical and "
//                     "analytical derivative wrt d"
//                  << it.first << " are differnt: \n    "
//                  << "numerical value = " << num_der
//                  << "\n    analytical value = " << der << ". \n";
//     if (std::abs(num_der - der) / std::abs(num_der) > 1e-6)
//       Moose::out << "checksoftMaxMinGradient the value of the numerical and "
//                     "analytical derivative wrt d"
//                  << it.first
//                  << " are differnt:\n    numerical value = " << num_der
//                  << "\n    analytical value = " << der << ". \n";
//   }
// }
