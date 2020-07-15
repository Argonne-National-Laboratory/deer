#pragma once

#include "nlFunBase.h"

class f_aL : public nlFunBase {

public:
  f_aL(const std::string &f_name, const Real &D);
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

protected:
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {"a"};
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp(0);
    return temp = {"eqedotc", "sVM"};
  };

private:
  const Real _n = 1. / 3.;
  Real calc_L(const io_maps_type &params) const;
  const Real _D;
};
