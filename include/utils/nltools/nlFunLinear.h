#pragma once

#include "nlFunBase.h"

class nlFunLinear : public nlFunBase {

public:
  nlFunLinear(const std::string &fname, const std::string &var_name,
              const Real m, const Real q);
  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override;
  io_maps_type
  computeVarGradient(const io_maps_type &x,
                     const io_maps_type &params,
                     const io_maps_type &x_old) const override;
  io_maps_type
  computeParamGradient(const io_maps_type &x,
                       const io_maps_type &params,
                       const io_maps_type &x_old) const override;

  void setMandQ(const Real &m, const Real &q) {
    _m = m;
    _q = q;
  };

protected:
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {};
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp = {};
    return temp;
  };
  const std::string _var_name;
  Real _m;
  Real _q;
};
