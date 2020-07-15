#pragma once

#include "equationTests.h"

equationTests::equationTests(const std::string &fname, const Real &a0,
                             const Real &b0, const Real &sigma_0,
                             const Real &S_thr, const Real &beta_exponent,
                             const Real &b_saturation, const Real &E_interface,
                             const bool &cavity_nucleation_on,
                             const bool &cavity_growth_on,
                             const bool &triaxial_cavity_growth_on)
    : nlFunBase(fname, requiredVar(), requiredParam()), _a0(a0), _b0(b0),
      _sigma_0(sigma_0), _S_thr(S_thr), _beta_exponent(beta_exponent),
      _b_saturation(b_saturation), _E_interface(E_interface),
      _cavity_nucleation_on(cavity_nucleation_on),
      _cavity_growth_on(cavity_growth_on),
      _triaxial_cavity_growth_on(triaxial_cavity_growth_on) {}

bool equationTests::nucleationIsActive(const io_maps_type &x,
                                       const io_maps_type &params,
                                       const io_maps_type &x_old) const {

  if (!_cavity_nucleation_on)
    return false;
  Real Tn_old = getValueFromMap(x_old, "Tn", "equationTests::computeValue.x");
  // if (Tn_old < 1)
  //   return false;

  Real uN_dot_old =
      getValueFromMap(params, "un_dot_old", "nucleationIsActive.params");
  Real uN_dot = getValueFromMap(params, "un_dot", "bdot::computeValue.params");
  //
  // if (uN_dot_old < 0 && uN_dot > 0)
  //   return false;

  Real b_old = getValueFromMap(x_old, "b", "bdot::computeValue.x_old");
  Real eqec = getValueFromMap(params, "eqec", "bdot::computeValue.params");
  Real uN = getValueFromMap(params, "un", "bdot::computeValue.params");

  Real Tn = getValueFromMap(x, "Tn", "bdot::computeValue.x");

  /* check saturation */
  if (b_old <= _b_saturation)
    // std::cout << "b_old <= _b_saturation" << std::endl;
    return false;

  /* nucleation has already been activated once*/
  if (getValueFromMap(params, "nucleation_above_threshold") == 1)
    return true;

  return false;

  // if (Tn_old <= 0)
  //   return false;
  //
  // if (Tn_old > 0 && uN_dot > 0)
  //   return true;

  // /*obvious compression state*/
  // if (Tn_old <= 0 && un_dot <= 0)
  //   return false;
  //
  // if (Tn_old <= 0 && !linearUpdate(x, params, x_old))
  //   return false;
  //
  // /*obvious tension state*/
  // if (Tn_old > 0 && uN_dot > 0)
  //   return true;
  //
  // /* if it is not obvious but un_dot is positive we need to assume
  // nucleation is
  //  * accitve*/
  // if (uN_dot > 0)
  //   return true;
  //
  // /*if we are here no nucleation*/
  // return false;
}

bool equationTests::nucleationAboveThreshold(const io_maps_type &params,
                                             const io_maps_type &x_old) const {

  if (!_cavity_nucleation_on)
    return false;

  Real Tn_old = getValueFromMap(x_old, "Tn", "equationTests::computeValue.x");
  Real eqec = getValueFromMap(params, "eqec", "bdot::computeValue.params");

  Real S = std::pow(Tn_old / _sigma_0, _beta_exponent) * eqec;
  if (S > _S_thr & Tn_old > 0)
    return true;

  return false;
}
bool equationTests::growthIsActive(const io_maps_type &x,
                                   const io_maps_type &params,
                                   const io_maps_type &x_old) const {

  if (!_cavity_growth_on)
    return false;

  Real uN_dot_old =
      getValueFromMap(params, "un_dot_old", "nucleationIsActive.params");
  Real uN_dot = getValueFromMap(params, "un_dot", "bdot::computeValue.params");

  // if (uN_dot_old < 0 && uN_dot > 0)
  //   return false;

  Real a_old = getValueFromMap(x_old, "a", "adot::growthIsActive");
  Real TN_old = getValueFromMap(x_old, "Tn", "adot::growthIsActive");
  Real uN = getValueFromMap(params, "un", "adot::growthIsActive");

  return true;
}

bool equationTests::triaxialGrowthIsActive(const io_maps_type &x,
                                           const io_maps_type &params,
                                           const io_maps_type &x_old) const {

  if (!_triaxial_cavity_growth_on)
    return false;

  Real TN_old = getValueFromMap(x_old, "Tn", "adot::triaxialGrowthIsActive");

  // if (TN_old < 0)
  //   return false;

  // Real sH = getValueFromMap(params, "sH", "adot::triaxialGrowthIsActive");
  // if (sH <= 0)
  //   return false;

  // Real uN = getValueFromMap(params, "un", "adot::triaxialGrowthIsActive");
  // Real uNOld =
  //     getValueFromMap(params, "un_old", "adot::triaxialGrowthIsActive");
  // if (std::copysign(1., uN) * std::copysign(1., uNOld) != 1.)
  //   return false;

  return true;
}

bool equationTests::linearUpdate(const io_maps_type &x,
                                 const io_maps_type &params,
                                 const io_maps_type &x_old) const {

  Real a = getValueFromMap(x, "a", "adot::growthIsActive");
  Real b = getValueFromMap(x, "b", "adot::growthIsActive");
  Real uN = getValueFromMap(params, "un", "bdot::computeValue.params");
  Real uN_old = getValueFromMap(params, "un_old", "bdot::computeValue.params");
  Real TN_old = getValueFromMap(x_old, "Tn", "bdot::computeValue.params");
  Real TNlinear =
      (uN - uN_old) * _E_interface * (1 - a / b) / (6 * _b0) + TN_old;

  return TNlinear > 0;
}
