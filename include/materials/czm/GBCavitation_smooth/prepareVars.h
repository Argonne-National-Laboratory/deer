#pragma once

#include "nlFunBase.h"

class prepareVars : public nlFunBase {

public:
  prepareVars(const std::string &f_name, const Real &a0);
  io_maps_type computeNewX(const io_maps_type &x,
                           const io_maps_type &x_old) const;
  io_maps_type computeNewXVarGradient(const io_maps_type &x,
                                      const io_maps_type &x_old) const;

  prepareVars::io_maps_type XNLFromComputedX(const io_maps_type &comp_x,
                                             const io_maps_type &x_old) const;

protected:
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {"a", "b"};
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp(0);
    return temp;
  };
  const Real _a0;
};
