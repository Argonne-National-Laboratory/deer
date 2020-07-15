#pragma once

#include "MathUtils.h"
#include "MooseError.h"
#include <unordered_map>

#define UNUSED(x) (void)(x)

class nlFunBase {

public:
  typedef std::unordered_map<std::string, Real> io_maps_type;
  typedef std::unordered_map<std::string, io_maps_type *> map_of_io_maps_type;
  typedef std::vector<std::string> str_vct_type;
  const io_maps_type _empty_map;
  nlFunBase(const std::string & /*f_name*/,
            const str_vct_type required_var = {},
            const str_vct_type required_param = {});
  nlFunBase(const std::string & /*f_name*/,
            const str_vct_type /*required_var = {}*/,
            const str_vct_type /*required_param = {}*/,
            const bool /*fun_need_old_x*/);

  virtual void prepare(const io_maps_type & /*x*/,
                       const io_maps_type & /*params*/,
                       const io_maps_type & /*x_old*/) {}

  bool nlFunNeedOldX() const { return _fun_need_old_x; }

  // std::string getFunName() const { return _f_name; };
  // Real getValue() const { return _f_value; };
  // io_maps_type getVarGradient() const { return _f_var_gradient; };
  // io_maps_type getParamGradient() const { return _f_param_gradient; };

  virtual Real computeValue(const io_maps_type &x, const io_maps_type &params,
                            const io_maps_type &x_old) const {
    UNUSED(x);
    UNUSED(x_old);
    UNUSED(params);
    notImplemented();
    return 0;
  }

  virtual io_maps_type computeVarGradient(const io_maps_type &x,
                                          const io_maps_type &params,
                                          const io_maps_type &x_old) const {
    UNUSED(x);
    UNUSED(x_old);
    UNUSED(params);
    notImplemented();
    return _empty_map;
  }

  virtual io_maps_type computeParamGradient(const io_maps_type &x,
                                            const io_maps_type &params,
                                            const io_maps_type &x_old) const {
    UNUSED(x);
    UNUSED(x_old);
    UNUSED(params);
    notImplemented();
    return _empty_map;
  }

  void computeAndSetValue(const io_maps_type &x, const io_maps_type &params,
                          const io_maps_type &x_old) {
    _f_value = computeValue(x, params, x_old);
  }

  void computeAndSetGradient(const io_maps_type &x, const io_maps_type &params,
                             const io_maps_type &x_old,
                             const bool &compute_var_garadient = true,
                             const bool &compute_param_garadient = false) {
    if (compute_var_garadient)
      _f_var_gradient = computeVarGradient(x, params, x_old);
    if (compute_param_garadient)
      _f_param_gradient = computeParamGradient(x, params, x_old);
  }

  void updateAndComputeAll(const io_maps_type &x, const io_maps_type &params,
                           const io_maps_type &x_old,
                           const bool &compute_var_garadient = true,
                           const bool &compute_param_garadient = false) {
    prepare(x, params, x_old);
    computeAndSetValue(x, params, x_old);
    computeAndSetGradient(x, params, x_old, compute_var_garadient,
                          compute_param_garadient);
  };

  io_maps_type computeNumericalGradient(const io_maps_type & /*x*/,
                                        const io_maps_type & /*params*/,
                                        const io_maps_type & /*x_old*/) const;

  void checkNumericalGradient(const io_maps_type & /*x*/,
                              const io_maps_type & /*params*/,
                              const io_maps_type &x_old) const;

  // // Overload + operator to add two nlFunBase objects.
  // nlFunBase operator+(const nlFunBase &b) const {
  //   std::string f_name;
  //   str_vct_type req_var, req_param;
  //   joinRequirements(b, "+", f_name, req_var, req_param);
  //   nlFunBase r(f_name, req_var, req_param);
  //
  //   Real va = this->_f_value;
  //   Real vb = b._f_value;
  //   r._f_value = va - vb;
  //   r._f_var_gradient =
  //       concatenateStringRealMap(this->_f_var_gradient, b._f_var_gradient);
  //
  //   for (const auto &it : r._f_var_gradient) {
  //     Real da_dx, db_dx;
  //     getPairedValues(this->_f_var_gradient, b._f_var_gradient, it.first,
  //     da_dx,
  //                     db_dx);
  //     setValueInMap(r._f_var_gradient, it.first, (da_dx + db_dx));
  //   }
  //
  //   r._f_param_gradient =
  //       concatenateStringRealMap(this->_f_param_gradient,
  //       b._f_param_gradient);
  //   for (const auto &it : r._f_param_gradient) {
  //     Real da_dx, db_dx;
  //     getPairedValues(this->_f_param_gradient, b._f_param_gradient, it.first,
  //                     da_dx, db_dx);
  //     setValueInMap(r._f_param_gradient, it.first, da_dx + db_dx);
  //   }
  //   return r;
  // }
  //
  // // Overload log .
  // nlFunBase logFun(const nlFunBase &a) const {
  //   std::string f_name = a.getFunName();
  //   io_maps_type f_v_grad = a.getVarGradient();
  //   io_maps_type f_p_grad = a.getParamGradient();
  //   nlFunBase r = a;
  //
  //   r._f_value = std::log(a.getValue());
  //   Real dlog_df = -1. / std::log(a.getValue());
  //
  //   r._f_var_gradient = chain(dlog_df, f_v_grad);
  //   r._f_param_gradient = chain(dlog_df, f_p_grad);
  //
  //   return r;
  // }
  //
  // nlFunBase operator-(const nlFunBase &b) const {
  //   std::string f_name;
  //   str_vct_type req_var, req_param;
  //   joinRequirements(b, "-", f_name, req_var, req_param);
  //   nlFunBase r(f_name, req_var, req_param);
  //
  //   Real va = this->_f_value;
  //   Real vb = b._f_value;
  //   r._f_value = va - vb;
  //   r._f_var_gradient =
  //       concatenateStringRealMap(this->_f_var_gradient, b._f_var_gradient);
  //
  //   for (const auto &it : r._f_var_gradient) {
  //     Real da_dx, db_dx;
  //     getPairedValues(this->_f_var_gradient, b._f_var_gradient, it.first,
  //     da_dx,
  //                     db_dx);
  //     setValueInMap(r._f_var_gradient, it.first, (da_dx - db_dx));
  //   }
  //
  //   r._f_param_gradient =
  //       concatenateStringRealMap(this->_f_param_gradient,
  //       b._f_param_gradient);
  //   for (const auto &it : r._f_param_gradient) {
  //     Real da_dx, db_dx;
  //     getPairedValues(this->_f_param_gradient, b._f_param_gradient, it.first,
  //                     da_dx, db_dx);
  //     setValueInMap(r._f_param_gradient, it.first, da_dx - db_dx);
  //   }
  //   return r;
  // }
  //
  // nlFunBase operator*(const nlFunBase &b) const {
  //   std::string f_name;
  //   str_vct_type req_var, req_param;
  //   joinRequirements(b, "*", f_name, req_var, req_param);
  //   nlFunBase r(f_name, req_var, req_param);
  //
  //   Real va = this->_f_value;
  //   Real vb = b._f_value;
  //   r._f_value = va * vb;
  //   r._f_var_gradient =
  //       concatenateStringRealMap(this->_f_var_gradient, b._f_var_gradient);
  //
  //   for (const auto &it : r._f_var_gradient) {
  //     Real da_dx, db_dx;
  //     getPairedValues(this->_f_var_gradient, b._f_var_gradient, it.first,
  //     da_dx,
  //                     db_dx);
  //     setValueInMap(r._f_var_gradient, it.first, (da_dx * vb + va * db_dx));
  //   }
  //
  //   r._f_param_gradient =
  //       concatenateStringRealMap(this->_f_param_gradient,
  //       b._f_param_gradient);
  //   for (const auto &it : r._f_param_gradient) {
  //     Real da_dx, db_dx;
  //     getPairedValues(this->_f_param_gradient, b._f_param_gradient, it.first,
  //                     da_dx, db_dx);
  //     setValueInMap(r._f_param_gradient, it.first, da_dx * vb + va * db_dx);
  //   }
  //   return r;
  // }
  //
  // nlFunBase operator/(const nlFunBase &b) const {
  //   std::string f_name;
  //   str_vct_type req_var, req_param;
  //   joinRequirements(b, "/", f_name, req_var, req_param);
  //   nlFunBase r(f_name, req_var, req_param);
  //
  //   Real va = this->_f_value;
  //   Real vb = b._f_value;
  //   Real vb2 = vb * vb;
  //   r._f_value = va / vb;
  //   r._f_var_gradient =
  //       concatenateStringRealMap(this->_f_var_gradient, b._f_var_gradient);
  //
  //   for (const auto &it : r._f_var_gradient) {
  //     Real da_dx, db_dx;
  //     getPairedValues(this->_f_var_gradient, b._f_var_gradient, it.first,
  //     da_dx,
  //                     db_dx);
  //     setValueInMap(r._f_var_gradient, it.first,
  //                   (da_dx * vb - db_dx * va) / vb2);
  //   }
  //
  //   r._f_param_gradient =
  //       concatenateStringRealMap(this->_f_param_gradient,
  //       b._f_param_gradient);
  //   for (const auto &it : r._f_param_gradient) {
  //     Real da_dx, db_dx;
  //     getPairedValues(this->_f_param_gradient, b._f_param_gradient, it.first,
  //                     da_dx, db_dx);
  //     setValueInMap(r._f_param_gradient, it.first,
  //                   (da_dx * vb - db_dx * va) / vb2);
  //   }
  //   return r;
  // }

protected:
  void notImplemented() const {
    mooseError(_f_name + ":: method not implemented ");
  };
  void getPairedValues(const io_maps_type &a, const io_maps_type &b,
                       const std::string &mapkey, Real &va, Real &vb) const {
    va = 0;
    vb = 0;
    if (checkMapKeyExist(a, mapkey))
      va = getValueFromMap(a, mapkey, "getPairedValues::a");
    if (checkMapKeyExist(b, mapkey))
      vb = getValueFromMap(b, mapkey, "getPairedValues::b");
  }

  io_maps_type chain(const Real &df_dx, const io_maps_type &dx_dy) const {
    io_maps_type df_dy;
    for (const auto &it : dx_dy)
      df_dy[it.first] = it.second * df_dx;

    return df_dy;
  }

  io_maps_type sumD(const io_maps_type &df_dx,
                    const io_maps_type &dg_dx) const {
    io_maps_type dfdx_plus_dgdx = concatenateStringRealMap(df_dx, dg_dx);

    for (const auto &it : dfdx_plus_dgdx) {
      Real v1, v2;
      getPairedValues(df_dx, dg_dx, it.first, v1, v2);
      setValueInMap(dfdx_plus_dgdx, it.first, v1 + v2);
    }
    return dfdx_plus_dgdx;
  }

  io_maps_type subtractD(const io_maps_type &df_dx,
                         const io_maps_type &dg_dx) const {
    io_maps_type dfdx_minus_dgdx = concatenateStringRealMap(df_dx, dg_dx);

    for (const auto &it : dfdx_minus_dgdx) {
      Real v1, v2;
      getPairedValues(df_dx, dg_dx, it.first, v1, v2);
      setValueInMap(dfdx_minus_dgdx, it.first, v1 - v2);
    }
    return dfdx_minus_dgdx;
  }

  io_maps_type f_divided_g_D(const Real &f, const io_maps_type &df_dx,
                             const Real &g, const io_maps_type &dg_dx) const {

    return chain(1. / (g * g), subtractD(chain(g, df_dx), chain(f, dg_dx)));
  }

  io_maps_type f_times_g_D(const Real &f, const io_maps_type &df_dx,
                           const Real &g, const io_maps_type &dg_dx) const {

    return sumD(chain(g, df_dx), chain(f, dg_dx));
  }

  io_maps_type f_power_n_D(const Real &f, const io_maps_type &df_dx,
                           const Real &n) const {

    return chain(n * std::pow(f, n - 1.), df_dx);
  }

  str_vct_type joinRequirementVector(const str_vct_type &v1,
                                     const str_vct_type &v2) const {
    str_vct_type req_vct = concatenateStringVector(v1, v2);
    std::sort(req_vct.begin(), req_vct.end());
    req_vct.erase(unique(req_vct.begin(), req_vct.end()), req_vct.end());
    return req_vct;
  }

  void joinRequirements(const nlFunBase &b, const std::string &s,
                        std::string f_name, str_vct_type req_var,
                        str_vct_type req_param) const {
    f_name = this->_f_name + s + b._f_name;
    req_var = joinRequirementVector(this->_required_var, b._required_var);
    req_param = joinRequirementVector(this->_required_param, b._required_param);
  }

  io_maps_type initMaps(const str_vct_type &var_names);

  void optimalDx(const Real & /*x*/, Real & /*x_plus_dx*/, Real & /*dx*/,
                 const Real &tol = 1e-6) const;

  Real numericDiff(const Real & /*f0*/, const Real & /*f0_plus_dx*/,
                   const Real & /*dx*/) const;
  Real getValueFromMap(const io_maps_type & /*map*/,
                       const std::string & /*key*/,
                       const std::string &map_name = "") const;
  bool checkMapKeyExist(const io_maps_type & /*map*/,
                        const std::string & /*key*/) const;

  void setValueInMap(io_maps_type & /*map*/, const std::string & /*key*/,
                     const Real & /*value*/,
                     const std::string &map_name = "") const;
  str_vct_type concatenateStringVector(const str_vct_type & /*va*/,
                                       const str_vct_type & /*vb*/) const;

  io_maps_type concatenateStringRealMap(const io_maps_type & /*map_a*/,
                                        const io_maps_type & /*map_b*/) const;

  const std::string _f_name;
  Real _f_value;
  io_maps_type _f_var_gradient;
  io_maps_type _f_param_gradient;
  const str_vct_type _required_var;
  const str_vct_type _required_param;
  const bool _fun_need_old_x;
  const Real _pi = libMesh::pi;

  void checkAllRequiredExists(const io_maps_type &map,
                              const str_vct_type &vec) const {
    for (const auto &required_key : vec)
      checkMapKeyExist(map, required_key);
  }

  void checkParamExist(const io_maps_type &map) const {
    checkAllRequiredExists(map, _required_param);
  };
  void checkVarExist(const io_maps_type &map) const {
    checkAllRequiredExists(map, _required_var);
  };
};

// class nlFunction : public nlFunBase {
//
// public:
//   nlFunction(const std::string &function_name,
//              const str_vct_type &nl_vars_names,
//              const str_vct_type &parameters_names,
//              const io_maps_type &nl_vars_current_values,
//              const io_maps_type &parameters_values);
//
//   std::string getFunctionName() const { return _function_name; };
//   str_vct_type getAllVarNames() const { return _nl_vars_names; };
//
//   Real getFunctionValue() const { return _function_value; };
//
//   io_maps_type getFunctionGradient() const { return _function_gradient;
//   };
//
//   Real getDfunDX(const std::string & /*dvar_name*/) const;
//
//   /// shortcuts to compute and set function's value and function's
//   gradient void computeAndSetFunctionValue(); void
//   computeAndSetFunctionGradient();
//
//
// protected:
//   const std::string _function_name;
//   const str_vct_type _nl_vars_names;
//   const str_vct_type _parameters_names;
//   const str_vct_type _vars_plus_params_names;
//   Real _function_value;
//   const io_maps_type &_nl_vars_current_values;
//   const io_maps_type &_parameters_values;
//   io_maps_type _vars_plus_param_values;
//   io_maps_type _function_gradient;
//
//   io_maps_type initMaps(const str_vct_type &var_names);
//
//   void optimalDx(const Real &x, Real &x_plus_dx, Real &dx);
//
//   Real numericDiff(const Real &f0, const Real &f0_plus_dx, const Real
//   &dx); Real getValueFromMap(const io_maps_type &map, const std::string
//   &key,
//                        const std::string &map_name = "");
//   bool checkMapKeyExist(const io_maps_type &map, const std::string &key);
//
//   void setValueInMap(io_maps_type &map, const std::string &key,
//                      const Real &value, const std::string &map_name =
//                      "");
//   str_vct_type concatenateStringVector(const str_vct_type &va,
//                                        const str_vct_type &vb);
//
//   io_maps_type concatenateStringRealMap(const io_maps_type &map_a,
//                                         const io_maps_type &map_b);
//
//   // These are the function that need to be sublcassed to define a
//   function and
//   // its gradient. These function should not be called unless doing
//   numerical
//   // differentiation
//   virtual Real
//   computeFunctionValue(const io_maps_type & /*_var_values*/,
//                        const io_maps_type & /*_parameters_values*/) const
//                        {
//     return 0;
//   };
//   virtual io_maps_type
//   computeFunctionGradient(const io_maps_type & /*_var_values*/,
//                           const io_maps_type & /*_parameters_values*/
//                           ) const {
//     io_maps_type a;
//     return a;
//   };
//
//   io_maps_type computeFunctionNumericalGradient() const;
//   Real getValFromMap(const io_maps_type &map, const std::string &key,
//                      const std::string &map_name = "") const;
//   void setValInMap(io_maps_type &map, const std::string &key, const Real
//   &val,
//                    const std::string &map_name = "") const;
//
//   io_maps_type initMaps(const str_vct_type &var_names) const;
// };
// class nlVariable : public nlFunction {
//
// public:
//   nlVariable(const std::string & /*var_name*/, const Real &
//   /*_old_value*/,
//              const str_vct_type & /*depenedent_var_names*/,
//              const str_vct_type & /*parameters_names*/,
//              const io_maps_type & /*nl_vars_current_values*/,
//              const io_maps_type & /*parameters_values*/,
//              const Real & /*_scale_factor*/, const Real & /*offset*/,
//              const io_maps_type & /*dcurrent_dnl_all_vars*/);
//
//   void updateComputedValue();
//   void updateComputedValueGradient();
//   void updateResidual();
//   void updateResidualGradient();
//
//   void setOffset(const Real & /*newoffset*/);
//   void setScaleFactor(const Real & /*new_scale_factor*/);
//
//   Real getNLVarValue() const { return _nl_var_value; }
//   void setNLVarValue(const Real & /*new_nl_var_value*/);
//   Real getResidual(const bool &scaled = true) const;
//   Real getResidualDerivative(const std::string &dvar_name,
//                              const bool &scaled = true) const;
//
//   /// return dx_dxnl = dscale&offset*_dmodified_dnl
//   Real getDvarDvarnl() const { return _dxcurr_dxnl; };
//   Real getDvarnlDvar() const { return _dxnl_dxcurr; };
//
// protected:
//   Real getDxcurrentDnl(const std::string & /*dvar_name*/) const;
//   const std::string _var_name;
//
//   /// the variable value returned by the NL solver
//   Real _nl_var_value = 0;
//   /// the current var value after scaling back
//   Real _current_value;
//   /// the old solution var value
//   Real _old_value;
//   /// derivative of teh current value wrt to the nl value and viceversa
//   Real _dxcurr_dxnl = 0;
//   Real _dxnl_dxcurr = 0;
//
//   /// the value computed by the function in real and nl units
//   Real _computed_value;
//   Real _computed_value_nl;
//   Real _dxcompnl_dxcomp = 0;
//   Real _dxcomp_dxcompnl = 0;
//
//   /// the scale factor and offset defined as _current = _nl * scale +
//   offset Real _scale_factor; Real _offset;
//
//   /// boolean flags, mostly helper
//   bool _has_scale_factor;
//   bool _has_offset;
//
//   /// a map containing the derivatives of all current variable wrt to the
//   nl
//   /// variable value
//   const io_maps_type &_dcurrent_dnl_all_vars;
//
//   Real _residual = 0;
//   Real _residual_nl = 0;
//   io_maps_type _residual_gradient;
//   io_maps_type _residual_gradient_nl;
//
//   bool hasScaleFactor() const;
//   bool hasOffset() const;
//
//   /// this function is applied to the nl variable before scaling to
//   compute
//   /// the currnet_value
//   virtual void nlVarModifyFun(Real & /*nl_var_value*/,
//                               Real & /*modified_nl_value*/, Real &
//                               /*dout_din*/, const bool &modified_to_nl =
//                               false) const;
//
//   void scaleAndOffset(Real & /*modified_value*/, Real & /*real_value*/,
//                       Real & /*din_dout*/,
//                       const bool &real_to_modified = false) const;
// };
