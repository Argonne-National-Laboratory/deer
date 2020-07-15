#pragma once

#include "vdot.h"

vdot::vdot(const std::string &fname, const v1dot &v1L, const v1dot &v1H,
           const v2dotLow &v2L, const v2dotHigh &v2H,
           const equationTests &eqTest, const int &max_type,
           const int &vdot_type, const int &triaxial_vdot_active,
           const Real &vdot_smooth_factor, const Real &theta_time_integration)
    : nlFunRateNoHistory(fname, theta_time_integration), _v1L(v1L), _v1H(v1H),
      _v2L(v2L), _v2H(v2H),
      _vmax("max_of_abs_signed", /*f = */ vdot_smooth_factor,
            /*max_of_abs = */ false,
            /*max_of_abs_signed = */ true,
            /*compute_min = */ false),
      _eqTest(eqTest), _pnorm_max("pnormMax", vdot_smooth_factor),
      _soft_max("softMax", vdot_smooth_factor), _max_type(max_type),
      _vdot_type(vdot_type), _triaxial_vdot_active(triaxial_vdot_active) {}

Real vdot::computeValue(const io_maps_type &x, const io_maps_type &params,
                        const io_maps_type &x_old) const {

  if (!_eqTest.growthIsActive(x, params, x_old))
    return 0;

  if (getValueFromMap(x, "Tn") == 0 && getValueFromMap(params, "eqedotc") == 0)
    return 0;
  io_maps_type f;

  Real v1l = 0, v1h = 0, v2l = 0, v2h = 0;

  if (_vdot_type == 1 || _vdot_type == 3) {
    v1l = _v1L.computeValue(x, params, x_old);
    if (_eqTest.triaxialGrowthIsActive(x, params, x_old)) {
      v2l = _v2L.computeValue(x, params, x_old);
    }
  }

  if (_vdot_type == 2 || _vdot_type == 3) {
    v1h = _v1H.computeValue(x, params, x_old);
    if (_eqTest.triaxialGrowthIsActive(x, params, x_old)) {
      v2h = _v2H.computeValue(x, params, x_old);
    }
  }

  Real VL = v1l + v2l;
  Real VH = v1h + v2h;

  if (std::isfinite(v1l) && std::isfinite(v2l) && std::isfinite(v1h) &&
      std::isfinite(v2h)) {

    Real s;
    Real max;
    if (std::abs(VL) >= std::abs(VH)) {
      s = std::copysign(1., VL);
      max = VL;
    } else {
      s = std::copysign(1., VH);
      max = VH;
    }

    if (_max_type == 0)
      return max;

    else if (_max_type == 1) /*soft_max*/ {
      return _soft_max.computeSoftMaxValue(std::abs(VL), std::abs(VH)) * s;
    } else if (_max_type == 2) /*hard*/ {
      f["vL"] = VL;
      f["vH"] = VH;
      return _vmax.computeValue(f, params, x_old);
    } else if (_max_type == 3) /*staggered*/ {
      Real vLdot_old = getValueFromMap(x_old, "VL_dot");
      Real vHdot_old = getValueFromMap(x_old, "VH_dot");
      if (std::abs(vLdot_old) >= std::abs(vHdot_old))
        return VL;
      else
        return VH;
    } else if (_max_type == 4) /*average*/ {
      return (VL + VH) * 0.5;
    } else {
      mooseError("don't klnow this type of max");
      return 0;
    }
  } else
    mooseError(_f_name + "value is not finite: V1L " + std::to_string(v1l) +
               " V2L " + std::to_string(v2l) + " V1H " + std::to_string(v1h) +
               " V2H " + std::to_string(v2h));
  return 0;
}

vdot::io_maps_type vdot::computeVarGradient(const io_maps_type &x,
                                            const io_maps_type &params,
                                            const io_maps_type &x_old) const {

  if (!_eqTest.growthIsActive(x, params, x_old))
    return _empty_map;

  io_maps_type f;
  map_of_io_maps_type f_grad;

  if (getValueFromMap(x, "Tn") == 0 && getValueFromMap(params, "eqedotc") == 0)
    return _empty_map;

  Real v1l = 0, v1h = 0, v2l = 0, v2h = 0;
  io_maps_type dv1L_dx = {}, dv2L_dx = {}, dv1H_dx = {}, dv2H_dx = {};
  if (_vdot_type == 1 || _vdot_type == 3) {
    v1l = _v1L.computeValue(x, params, x_old);
    dv1L_dx = _v1L.computeVarGradient(x, params, x_old);
    if (_eqTest.triaxialGrowthIsActive(x, params, x_old)) {
      v2l = _v2L.computeValue(x, params, x_old);
      dv2L_dx = _v2L.computeVarGradient(x, params, x_old);
    }
  }

  if (_vdot_type == 2 || _vdot_type == 3) {
    v1h = _v1H.computeValue(x, params, x_old);
    dv1H_dx = _v1H.computeVarGradient(x, params, x_old);
    if (_eqTest.triaxialGrowthIsActive(x, params, x_old)) {
      v2h = _v2H.computeValue(x, params, x_old);
      dv2H_dx = _v2H.computeVarGradient(x, params, x_old);
    }
  }

  Real VL = v1l + v2l;
  Real VH = v1h + v2h;

  io_maps_type dVL_dx = sumD(dv1L_dx, dv2L_dx);
  io_maps_type dVH_dx = sumD(dv1H_dx, dv2H_dx);
  io_maps_type dmax_dx;
  if (std::isfinite(v1l) && std::isfinite(v2l) && std::isfinite(v1h) &&
      std::isfinite(v2h)) {

    Real s;
    if (std::abs(VL) >= std::abs(VH)) {
      s = std::copysign(1., VL);
      dmax_dx = dVL_dx;
    } else {
      s = std::copysign(1., VH);
      dmax_dx = dVH_dx;
    }

    if (_max_type == 0)
      return dmax_dx;

    else if (_max_type == 1) {
      io_maps_type dvi_dfi =
          _soft_max.computeSoftMaxGradient(std::abs(VL), std::abs(VH));

      Real dmax_dvl = dvi_dfi.find("a")->second * std::copysign(1., VL);
      Real dmax_dvh = dvi_dfi.find("b")->second * std::copysign(1., VH);
      dVL_dx = chain(dmax_dvl, dVL_dx);
      dVH_dx = chain(dmax_dvh, dVH_dx);
      io_maps_type dvdx = sumD(dVL_dx, dVH_dx);
      dvdx = chain(s, dvdx);

      return dvdx;
    } else if (_max_type == 2) {
      f["vL"] = VL;
      f["vH"] = VH;
      map_of_io_maps_type dVi_dxi;
      dVi_dxi["vL"] = &dVL_dx;
      dVi_dxi["vH"] = &dVH_dx;
      io_maps_type dV_dxi = _vmax.computeVarGradient(f, dVi_dxi, params);

      return dV_dxi;

    } else if (_max_type == 3) {
      Real vLdot_old = getValueFromMap(x_old, "VL_dot");
      Real vHdot_old = getValueFromMap(x_old, "VH_dot");
      if (std::abs(vLdot_old) >= std::abs(vHdot_old))
        return dVL_dx;
      else
        return dVH_dx;
    } else if (_max_type == 4) {
      return chain(0.5, sumD(dVL_dx, dVH_dx));
    } else {
      mooseError("don't klnow this type of max");
      return _empty_map;
    }

  } else {
    mooseError("v1l or v2l are not finite");
    return _empty_map;
  }
}

vdot::io_maps_type vdot::computeParamGradient(const io_maps_type &x,
                                              const io_maps_type &params,
                                              const io_maps_type &x_old) const {
  if (!_eqTest.nucleationIsActive(x, params, x_old) &&
      !_eqTest.growthIsActive(x, params, x_old))
    return _empty_map;

  io_maps_type f;
  map_of_io_maps_type f_grad;
  return _empty_map;
}
