#pragma once

#include "TdotBase.h"

TdotBase::TdotBase(const std::string &fname,
                   const std::string &displacement_name,
                   const std::string &traction_name, const Real &E,
                   const Real &interface_thickness,
                   const Real &theta_time_integration)

    : nlFunRate(fname, traction_name, theta_time_integration),
      _displacement_name(displacement_name),
      _displacement_dot_name(displacement_name + "_dot"),
      _traction_name(traction_name), _E(E),
      _interface_thickness(interface_thickness) {}

Real TdotBase::computeValue(const io_maps_type &x, const io_maps_type &params,
                            const io_maps_type &x_old) const {

  Real ux_dot = getValueFromMap(params, _displacement_dot_name,
                                "TdotBase::computeValue.x");

  Real T_dot = (ux_dot - getSecondTermValue(x, params, x_old)) *
               Stiffness(x, params, x_old);

  return T_dot;
}

TdotBase::io_maps_type
TdotBase::computeVarGradient(const io_maps_type &x, const io_maps_type &params,
                             const io_maps_type &x_old) const {

  io_maps_type f_grad;
  Real ux_dot = getValueFromMap(params, _displacement_dot_name,
                                "TdotBase::computeValue.x");

  Real C = Stiffness(x, params, x_old);
  Real dTdot_dC = ux_dot - getSecondTermValue(x, params, x_old);
  io_maps_type dC_dx = StiffnessVarGradient(x, params, x_old);
  f_grad = chain(dTdot_dC, dC_dx);

  Real dTdot_dsecond_term = -C;
  io_maps_type dsecond_term_dx = getSecondTermVarGradient(x, params, x_old);
  dsecond_term_dx = chain(dTdot_dsecond_term, dsecond_term_dx);

  f_grad = sumD(f_grad, dsecond_term_dx);

  return f_grad;
}

TdotBase::io_maps_type
TdotBase::computeParamGradient(const io_maps_type &x,
                               const io_maps_type &params,
                               const io_maps_type &x_old) const {
  io_maps_type f_grad;
  Real ux_dot = getValueFromMap(params, _displacement_dot_name,
                                "TdotBase::computeValue.x");

  Real C = Stiffness(x, params, x_old);
  Real dTdot_dC = ux_dot - getSecondTermValue(x, params, x_old);
  io_maps_type dC_dp = StiffnessParamGradient(x, params, x_old);
  f_grad = chain(dTdot_dC, dC_dp);

  Real dTdot_dsecond_term = -C;
  io_maps_type dsecond_term_dp = getSecondTermParamGradient(x, params, x_old);
  dsecond_term_dp = chain(dTdot_dsecond_term, dsecond_term_dp);

  f_grad = sumD(f_grad, dsecond_term_dp);

  io_maps_type dTdot_dudot;
  dTdot_dudot[_displacement_dot_name] = C;

  f_grad = sumD(f_grad, dTdot_dudot);

  return f_grad;
}

Real TdotBase::Stiffness(const io_maps_type &x, const io_maps_type &params,
                         const io_maps_type &x_old) const {
  UNUSED(params);
  UNUSED(x_old);
  Real a = getValueFromMap(x, "a", "TdotBase::computeValue.x");
  Real b = getValueFromMap(x, "b", "TdotBase::computeValue.x");

  return _E * (1 - a / b) / (_interface_thickness);
}

TdotBase::io_maps_type
TdotBase::StiffnessVarGradient(const io_maps_type &x,
                               const io_maps_type &params,
                               const io_maps_type &x_old) const {
  UNUSED(params);
  UNUSED(x_old);
  Real a = getValueFromMap(x, "a", "TdotBase::computeValue.x");
  Real b = getValueFromMap(x, "b", "TdotBase::computeValue.x");
  io_maps_type dCdx;
  dCdx["a"] = -_E / (b * _interface_thickness);
  dCdx["b"] = _E * a / (b * b * _interface_thickness);
  return dCdx;
}
