#pragma once

#include "nlFunBase.h"

class softMax : public nlFunBase {

public:
  softMax(const std::string & /*f_name*/, const Real &f = 50);

  // virtual void prepare(const io_maps_type &x,
  //                      const io_maps_type &params) override {
  //   updateX2useAndParam(x);
  // };

  Real computeSoftMaxValue(const Real & /*a*/, const Real & /*b*/) const;

  io_maps_type computeSoftMaxGradient(const Real & /*a*/,
                                      const Real & /*b*/) const;

protected:
  const Real _f;

  virtual str_vct_type requiredVar() {
    str_vct_type temp(0);
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp(0);
    return temp;
  };
};
