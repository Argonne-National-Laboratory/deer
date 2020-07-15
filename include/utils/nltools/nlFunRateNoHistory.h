#pragma once

#include "nlFunBase.h"

class nlFunRateNoHistory : public nlFunBase {

public:
  nlFunRateNoHistory(const std::string &fname,
                     const Real &theta_time_integration);

  virtual Real computeTimeIntegral(const io_maps_type &x,
                                   const io_maps_type &params,
                                   const io_maps_type &x_old,
                                   const Real &dt) const;

  virtual io_maps_type computeTimeIntegralVarGradient(
      const io_maps_type &x, const io_maps_type &params,
      const io_maps_type &x_old, const Real &dt) const;

  virtual io_maps_type computeTimeIntegralParamGradient(
      const io_maps_type &x, const io_maps_type &params,
      const io_maps_type &x_old, const Real &dt) const;

protected:
  const std::string _var_name;
  const Real _theta_implicit;
  const Real _theta_explicit;
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {};
    // mooseError("nlFunRateNoHistory::requiredVar has not been overridden,
    // something is
    // "
    //            "wrong if you need to call this method!!!!");
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp = {};
    // mooseError("nlFunRateNoHistory::requiredParam has not been overridden,
    // something is "
    //            "if you need to call this method!!!!");
    return temp;
  };
};
