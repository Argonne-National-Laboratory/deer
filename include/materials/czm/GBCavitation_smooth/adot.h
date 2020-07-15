#pragma once

#include "equationTests.h"
#include "nlFunRate.h"
#include "vdot.h"

class adot : public nlFunRate {

public:
  adot(const std::string &f_name, const vdot &vdotFun, const Real &h,
       const Real &a0, const equationTests &eqTest,
       const Real &theta_time_integration);

  // adot(const std::string &f_name, const Real &h, const Real &n, const Real
  // &D);

  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override;
  io_maps_type computeVarGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const override;
  io_maps_type computeParamGradient(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old) const override;

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

  Real dfdaFun(const io_maps_type &x, const io_maps_type &params,
               const io_maps_type &x_old) const;

  Real dfdvdotFun(const io_maps_type &x) const;

  const vdot &_vdot;
  const Real _h;
  const Real _a0;
  const equationTests &_eqTest;
};
