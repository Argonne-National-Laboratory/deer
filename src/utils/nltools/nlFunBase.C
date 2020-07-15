
#include "nlFunBase.h"

nlFunBase::nlFunBase(const std::string &f_name,
                     const str_vct_type required_var /*= {}*/,
                     const str_vct_type required_param /*= {}*/)
    : _f_name(f_name), _required_var(required_var),
      _required_param(required_param), _fun_need_old_x(false) {}

nlFunBase::nlFunBase(const std::string &f_name,
                     const str_vct_type required_var /*= {}*/,
                     const str_vct_type required_param /*= {}*/,
                     const bool fun_need_old_x)
    : _f_name(f_name), _required_var(required_var),
      _required_param(required_param), _fun_need_old_x(fun_need_old_x) {}

nlFunBase::io_maps_type nlFunBase::initMaps(const str_vct_type &var_names) {
  io_maps_type map;
  for (auto &name : var_names)
    map[name] = 0;
  return map;
}

void nlFunBase::optimalDx(const Real &x, Real &x_plus_dx, Real &dx,
                          const Real &tol) const {
  if (x != 0)
    dx = x * tol;
  else
    dx = tol;
  x_plus_dx = x + dx;
}

Real nlFunBase::numericDiff(const Real &f0, const Real &f0_plus_dx,
                            const Real &dx) const {
  return (f0_plus_dx - f0) / dx;
}

Real nlFunBase::getValueFromMap(const io_maps_type &map, const std::string &key,
                                const std::string &map_name) const {
  auto it = map.find(key);
  if (it != map.end())
    return it->second;
  else
    mooseError("getValueFromMap:: can't fine the key " + key +
               " in the map named " + map_name);
}
bool nlFunBase::checkMapKeyExist(const io_maps_type &map,
                                 const std::string &key) const {
  auto it = map.find(key);
  if (it != map.end())
    return true;
  else
    return false;
}

void nlFunBase::setValueInMap(io_maps_type &map, const std::string &key,
                              const Real &value,
                              const std::string &map_name) const {
  auto it = map.find(key);
  if (it != map.end())
    it->second = value;
  else
    mooseError("setValueInMap:: can't fine the key " + key +
               " in the map named " + map_name);
}

nlFunBase::str_vct_type
nlFunBase::concatenateStringVector(const str_vct_type &va,
                                   const str_vct_type &vb) const {
  str_vct_type vout;
  vout.reserve(va.size() + vb.size()); // preallocate memory
  vout.insert(vout.end(), va.begin(), va.end());
  vout.insert(vout.end(), vb.begin(), vb.end());
  return vout;
}

nlFunBase::io_maps_type
nlFunBase::concatenateStringRealMap(const io_maps_type &map_a,
                                    const io_maps_type &map_b) const {
  io_maps_type map_out;
  map_out.clear();
  map_out.insert(map_a.begin(), map_a.end());
  map_out.insert(map_b.begin(), map_b.end());
  return map_out;
}

nlFunBase::io_maps_type
nlFunBase::computeNumericalGradient(const io_maps_type &x,
                                    const io_maps_type &parameters,
                                    const io_maps_type &x_old) const {

  const Real fval0 = computeValue(x, parameters, x_old);
  Real tol = 1e-6;
  Real dx = 0;
  Real dp = 0;
  Real fval_dx = 0;
  Real fval_dp = 0;
  Real df_dx = 0;
  Real df_dp = 0;

  io_maps_type function_numerical_gradinet;
  function_numerical_gradinet.clear();
  for (const auto &it : x) {
    io_maps_type x_temp = x;
    Real x_curr = getValueFromMap(x_temp, it.first, "x_temp");
    Real new_x_curr = 0;

    optimalDx(x_curr, new_x_curr, dx, tol);
    setValueInMap(x_temp, it.first, new_x_curr, "x_temp");
    fval_dx = computeValue(x_temp, parameters, x_old);
    df_dx = numericDiff(fval0, fval_dx, dx);
    function_numerical_gradinet[it.first] = df_dx;
  }
  for (auto &it : parameters) {
    io_maps_type p_temp = parameters;
    Real p_curr = getValueFromMap(p_temp, it.first, "p_temp");

    Real new_p_curr = 0;
    optimalDx(p_curr, new_p_curr, dp, tol);
    setValueInMap(p_temp, it.first, new_p_curr, "p_temp");
    fval_dp = computeValue(x, p_temp, x_old);
    df_dp = numericDiff(fval0, fval_dp, dp);
    function_numerical_gradinet[it.first] = df_dp;
  }

  return function_numerical_gradinet;
}

void nlFunBase::checkNumericalGradient(const io_maps_type &x,
                                       const io_maps_type &parameters,
                                       const io_maps_type &x_old) const {

  const Real rel_tol = 1e-5;

  io_maps_type function_numerical_gradinet =
      computeNumericalGradient(x, parameters, x_old);
  io_maps_type function_var_gradinet = computeVarGradient(x, parameters, x_old);
  io_maps_type function_param_gradinet =
      computeParamGradient(x, parameters, x_old);
  io_maps_type function_gradinet =
      concatenateStringRealMap(function_var_gradinet, function_param_gradinet);
  for (const auto &it : x) {

    Real num_der = getValueFromMap(function_numerical_gradinet, it.first,
                                   "function_numerical_gradinet");
    Real der = 0;
    if (checkMapKeyExist(function_var_gradinet, it.first))
      der =
          getValueFromMap(function_var_gradinet, it.first, "function_gradinet");

    Real ref = std::abs(num_der);
    if (num_der == 0)
      ref = 1;

    Moose::out << "VARIABLE " << it.first << "**********\n";
    if (std::log10(std::abs(num_der - der) / ref) > std::log10(rel_tol))
      Moose::out
          << "FAIL!!!\n the numerical and analytical derivative of function "
          << _f_name << " wrt " << it.first
          << " are differnt:\n    numerical value = " << num_der
          << "\n    analytical value = " << der << " rel difference "
          << std::abs(num_der - der) / ref << ". \n";
    else
      Moose::out
          << "PASS!!!\n   the numerical and analytical derivative of function "
          << _f_name << " wrt " << it.first
          << " are the same :\n    numerical value = " << num_der
          << "\n    analytical value = " << der << " rel difference "
          << std::abs(num_der - der) / ref << ". \n";
  }
  for (const auto &it : parameters) {

    Real num_der = getValueFromMap(function_numerical_gradinet, it.first,
                                   "function_numerical_gradinet");
    Real der = 0;
    if (checkMapKeyExist(function_param_gradinet, it.first))
      der = getValueFromMap(function_param_gradinet, it.first,
                            "function_gradinet");
    Real ref = std::abs(num_der);
    if (num_der == 0)
      ref = 1;

    Moose::out << "PARAMETER " << it.first << "**********\n";
    if (std::log10(std::abs(num_der - der) / ref) > std::log10(rel_tol))
      Moose::out
          << "FAIL!!!\n the numerical and analytical derivative of function "
          << _f_name << " wrt " << it.first
          << " are differnt:\n    numerical value = " << num_der
          << "\n    analytical value = " << der << " rel difference "
          << std::abs(num_der - der) / ref << ". \n";
    else
      Moose::out
          << "PASS!!!\n   the numerical and analytical derivative of function "
          << _f_name << " wrt " << it.first
          << " are the same :\n    numerical value = " << num_der
          << "\n    analytical value = " << der << " rel difference "
          << std::abs(num_der - der) / ref << ". \n";
  }
}

// nlFunction::nlFunction(const std::string &function_name,
//                        const str_vct_type &nl_vars_names,
//                        const str_vct_type &parameters_names,
//                        const io_maps_type &nl_vars_current_values,
//                        const io_maps_type &parameters_values)
//     : nlFunBase(), _function_name(function_name),
//     _nl_vars_names(nl_vars_names),
//       _parameters_names(parameters_names),
//       _vars_plus_params_names(
//           concatenateStringVector(nl_vars_names,
//           parameters_names)),
//       _function_value(0), _nl_vars_current_values(nl_vars_current_values),
//       _parameters_values(parameters_values),
//       _vars_plus_param_values(initMaps(_vars_plus_params_names)),
//       _function_gradient(initMaps(_vars_plus_params_names)) {}
//
// void nlFunction::computeAndSetFunctionValue() {
//   _function_value =
//       computeFunctionValue(_nl_vars_current_values, _parameters_values);
// };
//
// void nlFunction::computeAndSetFunctionGradient() {
//   _function_gradient =
//       computeFunctionGradient(_nl_vars_current_values, _parameters_values);
// };
//
// Real nlFunction::getDfunDX(const std::string &dvar_name) const {
//   return getValueFromMap(_function_gradient, dvar_name,
//                                     "_function_gradient");
// }
//
// Real nlFunction::getValFromMap(const io_maps_type &map, const std::string
// &key,
//                                const std::string &map_name) const {
//   return getValueFromMap(map, key, map_name);
// }
//
// void nlFunction::setValInMap(io_maps_type &map, const std::string &key,
//                              const Real &val,
//                              const std::string &map_name) const {
//   return setValueInMap(map, key, val, map_name);
// }
//
// io_maps_type nlFunction::initMaps(const str_vct_type &var_names) const {
//   return initMaps(var_names);
// };
//
// nlVariable::nlVariable(const std::string &var_name, const Real &old_value,
//                        const str_vct_type &depenedent_var_names,
//                        const str_vct_type &parameters_names,
//                        const io_maps_type &nl_vars_current_values,
//                        const io_maps_type &parameters_values,
//                        const Real &scale_factor, const Real &offset,
//                        const io_maps_type &dcurrent_dnl_all_vars)
//
//     : nlFunction("fun_" + var_name, depenedent_var_names, parameters_names,
//                  nl_vars_current_values, parameters_values),
//       _var_name(var_name), _current_value(0), _old_value(old_value),
//       _computed_value(0), _computed_value_nl(0), _scale_factor(scale_factor),
//       _offset(offset), _has_scale_factor(hasScaleFactor()),
//       _has_offset(hasOffset()),
//       _dcurrent_dnl_all_vars(dcurrent_dnl_all_vars),
//       _residual_gradient(initMaps(depenedent_var_names)),
//       _residual_gradient_nl(initMaps(depenedent_var_names)) {}
//
// void nlVariable::setScaleFactor(const Real &new_scale_factor) {
//   _scale_factor = new_scale_factor;
//   _has_scale_factor = hasScaleFactor();
// }
//
// Real nlVariable::getDxcurrentDnl(const std::string &dvar_name) const {
//   return getValueFromMap(_dcurrent_dnl_all_vars, dvar_name,
//                                     "_dcurrent_dnl_all_vars");
// }
//
// void nlVariable::setOffset(const Real &newoffset) {
//   _offset = newoffset;
//   _has_offset = hasOffset();
// }
//
// void nlVariable::updateComputedValue() {
//   computeAndSetFunctionValue();
//   _computed_value = getFunctionValue();
//   scaleAndOffset(_computed_value_nl, _computed_value, _dxcompnl_dxcomp,
//                  /*real_to_nl = */ true);
//   _dxcomp_dxcompnl = 1 / _dxcompnl_dxcomp;
// }
//
// void nlVariable::updateResidual() {
//   _residual = _current_value - _computed_value;
//   _residual_nl = _nl_var_value - _computed_value_nl;
// }
//
// Real nlVariable::getResidual(const bool &scaled) const {
//   if (!scaled)
//     return _residual;
//   else
//     return _residual_nl;
// }
//
// void nlVariable::updateResidualGradient() {
//   computeAndSetFunctionGradient();
//   for (auto &name : _nl_vars_names) {
//     Real dR_dxcurrentj = getDfunDX(name);
//     Real dxcurrentj_dxnlj = getDxcurrentDnl(name);
//     Real dRnl_dxnlj = _dxcompnl_dxcomp * dR_dxcurrentj * dxcurrentj_dxnlj;
//     if (name == _var_name) {
//       dR_dxcurrentj += 1;
//       dRnl_dxnlj += 1;
//     }
//     setValueInMap(_residual_gradient, name, dR_dxcurrentj,
//                              "_residual_gradient");
//     setValueInMap(_residual_gradient_nl, name, dRnl_dxnlj,
//                              "_residual_gradient_nl");
//   }
// }
//
// Real nlVariable::getResidualDerivative(const std::string &dvar_name,
//                                        const bool &scaled) const {
//   if (!scaled)
//     return getValueFromMap(_residual_gradient, dvar_name,
//                                       "_residual_gradient");
//   else
//     return getValueFromMap(_residual_gradient_nl, dvar_name,
//                                       "_residual_gradient_nl");
// }
//
// void nlVariable::nlVarModifyFun(Real &nl_var_value, Real &modified_nl_value,
//                                 Real &dout_din,
//                                 const bool &modified_to_nl) const {
//
//   if (!modified_to_nl) {
//     // modified_nl_value is some function of nl_var_value
//     modified_nl_value = nl_var_value;
//     // the derivative of modified_nl_value wrt nl_var_value
//     dout_din = 1;
//   } else {
//     //  nl_var_value is some function of modified_nl_value
//     nl_var_value = modified_nl_value;
//     dout_din = 1;
//   }
// }
//
// void nlVariable::scaleAndOffset(Real &modified_value, Real &real_value,
//                                 Real &dout_din,
//                                 const bool &real_to_modified) const {
//   if (!real_to_modified) { /// nl to real
//     // some function of modified_value
//     real_value = modified_value * _scale_factor + _offset;
//     // the derivative of real_value wrt modified_value
//     dout_din = _scale_factor;
//   } else { // some function of real_value
//     modified_value = (real_value - _offset) / _scale_factor;
//     // the derivative of modified_value wrt real_value
//     dout_din = 1 / _scale_factor;
//   }
// }
//
// void nlVariable::setNLVarValue(const Real &new_nl_var_value) {
//   _nl_var_value = new_nl_var_value;
//   Real mx, dmx_dxnl, dxnl_dmx;
//   Real dxcurr_dmx, dmx_dxcurr;
//   Real temp_trash;
//   nlVarModifyFun(_nl_var_value, mx, dmx_dxnl);
//   scaleAndOffset(mx, _current_value, dxcurr_dmx);
//   _dxcurr_dxnl = dxcurr_dmx * dmx_dxnl;
//   _dxnl_dxcurr = 1 / _dxcurr_dxnl;
// }
//
// bool nlVariable::hasScaleFactor() const {
//   if (_scale_factor != 1)
//     return true;
//   else
//     return false;
// }
//
// bool nlVariable::hasOffset() const {
//   if (_offset != 0)
//     return true;
//   else
//     return false;
// }
