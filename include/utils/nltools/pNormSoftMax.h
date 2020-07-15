#pragma once

#include "nlFunBase.h"

class pNormSoftMax : public nlFunBase {

public:
  pNormSoftMax(const std::string & /*f_name*/, const Real &p = 30);

  // virtual void prepare(const io_maps_type &x,
  //                      const io_maps_type &params) override {
  //   updateX2useAndParam(x);
  // };

  Real computeSoftMaxValue(const io_maps_type &x) const;

  io_maps_type computeSoftMaxGradient(const io_maps_type &x) const;

protected:
  const Real _p;
  virtual str_vct_type requiredVar() {
    str_vct_type temp(0);
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp(0);
    return temp;
  };
};
