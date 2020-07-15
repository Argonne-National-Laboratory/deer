#pragma once

#include "nlFunBase.h"

class logistic : public nlFunBase {

public:
  logistic(const std::string & /*f_name*/, const Real &L = 1, const Real &k = 1,
           const Real &x0 = 0);

  // virtual void prepare(const io_maps_type &x,
  //                      const io_maps_type &params) override {
  //   updateX2useAndParam(x);
  // };

  Real computeLogisticValue(const Real & /*x*/) const;

  io_maps_type computeLogisticGradient(const Real & /*x*/) const;

protected:
  const Real _L;
  const Real _k;
  const Real _x0;

  virtual str_vct_type requiredVar() {
    str_vct_type temp(0);
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp(0);
    return temp;
  };
};
