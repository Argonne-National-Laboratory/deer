#pragma once

#include "nlFunBase.h"

class nlFunConstant : public nlFunBase {

public:
  nlFunConstant(const std::string & /*f_name*/, const Real & /*C*/);
  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override;
  io_maps_type
  computeVarGradient(const io_maps_type &x, const io_maps_type &params,
                     const io_maps_type &x_old) const override;
  io_maps_type
  computeParamGradient(const io_maps_type &x,
                       const io_maps_type &params,
                       const io_maps_type &x_old) const override;

protected:
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {};
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp = {};
    return temp;
  };

  Real _C = 0;
};
