#pragma once

#include "TSdot.h"

TSdot::TSdot(const std::string &fname, const std::string &displacement_name,
             const std::string &traction_name, const Real &E,
             const Real &interface_thickness, const Real &eta,
             const Real &theta_time_integration)

    : TdotBase(fname, displacement_name, traction_name, E, interface_thickness,
               theta_time_integration),
      _eta(eta) {}

Real TSdot::getSecondTermValue(const io_maps_type &x,
                               const io_maps_type &params,
                               const io_maps_type &x_old) const {
  UNUSED(params);
  UNUSED(x_old);
  Real Ts = getValueFromMap(x, _traction_name, "TSdot::getSecondTermValue");
  return Ts / _eta;
}

TSdot::io_maps_type
TSdot::getSecondTermVarGradient(const io_maps_type &x,
                                const io_maps_type &params,
                                const io_maps_type &x_old) const {
  UNUSED(x);
  UNUSED(params);
  UNUSED(x_old);

  io_maps_type dsec_dx;
  dsec_dx[_traction_name] = 1 / _eta;
  return dsec_dx;
}
TSdot::io_maps_type
TSdot::getSecondTermParamGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const {

  UNUSED(x);
  UNUSED(params);
  UNUSED(x_old);
  return _empty_map;
}
