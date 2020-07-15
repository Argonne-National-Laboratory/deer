#pragma once

#include "nlFunBase.h"

class scaleFactors : public nlFunBase {

public:
  scaleFactors(const std::string &f_name);
  io_maps_type computeVarScaleFactor(const io_maps_type &x_old) const;
  io_maps_type computeVarScaleFactorGradient(const io_maps_type &x_old) const {
    UNUSED(x_old);
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
