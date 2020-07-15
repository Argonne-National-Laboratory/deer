#pragma once

#include "TNdot.h"

TNdot::TNdot(const std::string &fname, const std::string &displacement_name,
             const std::string &traction_name, const Real &E,
             const Real &interface_thickness, const vdot &vDot,
             const Real &E_penalty, const Real &theta_time_integration)
    : TdotBase(fname, displacement_name, traction_name, E, interface_thickness,
               theta_time_integration),
      _vdot(vDot), _E_Penalty(E_penalty) {}

Real TNdot::Stiffness(const io_maps_type &x, const io_maps_type &params,
                      const io_maps_type &x_old) const {
  Real C = TdotBase::Stiffness(x, params, x_old);
  Real un = getValueFromMap(params, _displacement_name, "TNdot::Stiffness");
  if (innerPenetration(un))
    C *= _E_Penalty;
  return C;
}

TNdot::io_maps_type
TNdot::StiffnessVarGradient(const io_maps_type &x, const io_maps_type &params,
                            const io_maps_type &x_old) const {

  io_maps_type dCdx = TdotBase::StiffnessVarGradient(x, params, x_old);
  Real un = getValueFromMap(params, _displacement_name, "TNdot::Stiffness");
  if (innerPenetration(un))
    dCdx = chain(_E_Penalty, dCdx);
  return dCdx;
}

bool TNdot::innerPenetration(const Real &un) const {
  return (un < -_interface_thickness);
}

Real TNdot::getSecondTermValue(const io_maps_type &x,
                               const io_maps_type &params,
                               const io_maps_type &x_old) const {
  UNUSED(params);
  UNUSED(x_old);
  Real b = getValueFromMap(x, "b", "TNdot::getSecondTermValue");
  return _vdot.computeValue(x, params, x_old) / (_pi * b * b);
}

TNdot::io_maps_type
TNdot::getSecondTermVarGradient(const io_maps_type &x,
                                const io_maps_type &params,
                                const io_maps_type &x_old) const {
  UNUSED(params);
  UNUSED(x_old);
  Real b = getValueFromMap(x, "b", "TNdot::getSecondTermValue");

  io_maps_type dsec_db;
  dsec_db["b"] =
      -2. * _vdot.computeValue(x, params, x_old) / (_pi * std::pow(b, 3.));
  Real dsec_dvdot = 1. / (_pi * b * b);
  io_maps_type dsec_dvdot_dvdot_dx =
      chain(dsec_dvdot, _vdot.computeVarGradient(x, params, x_old));

  io_maps_type d_second_term_dx = sumD(dsec_db, dsec_dvdot_dvdot_dx);
  return d_second_term_dx;
}
TNdot::io_maps_type
TNdot::getSecondTermParamGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const {
  UNUSED(params);
  UNUSED(x_old);
  Real b = getValueFromMap(x, "b", "TNdot::getSecondTermValue");
  Real dsec_dvdot = 1. / (_pi * b * b);
  return chain(dsec_dvdot, _vdot.computeParamGradient(x, params, x_old));
}
