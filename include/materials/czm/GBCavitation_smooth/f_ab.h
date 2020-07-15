#pragma once

#include "nlFunBase.h"

class f_ab : public nlFunBase {

public:
  f_ab(const std::string &f_name);
  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override;
  io_maps_type
  computeVarGradient(const io_maps_type &x,
                     const io_maps_type &params,
                     const io_maps_type &x_old) const override;
  io_maps_type
  computeParamGradient(const io_maps_type &x,
                       const io_maps_type &params,
                       const io_maps_type &x_old) const override {
    UNUSED(x);
    UNUSED(x_old);
    UNUSED(params);
    return _empty_map;
  };

protected:
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {"a", "b"};
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp(0);
    return temp;
  };
};
