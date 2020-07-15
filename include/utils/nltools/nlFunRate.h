#pragma once

#include "nlFunRateNoHistory.h"

class nlFunRate : public nlFunRateNoHistory {

public:
  nlFunRate(const std::string &fname, const std::string &var_name,
            const Real &theta_time_integration);

  virtual Real computeTimeIntegral(const io_maps_type &x,
                                   const io_maps_type &params,
                                   const io_maps_type &x_old,
                                   const Real &dt) const;

protected:
  const std::string _var_name;
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {};
    // mooseError("nlFunRate::requiredVar has not been overridden, something is
    // "
    //            "wrong if you need to call this method!!!!");
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp = {};
    // mooseError("nlFunRate::requiredParam has not been overridden, something
    // is "
    //            "if you need to call this method!!!!");
    return temp;
  };
};
