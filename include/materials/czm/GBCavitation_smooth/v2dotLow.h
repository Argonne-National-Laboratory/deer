#pragma once

#include "v2dotBase.h"

class v2dotLow : public v2dotBase {

public:
  v2dotLow(const std::string &fname, const Real &n, const Real &h,
           const Real &theta_time_integration, const equationTests &eqTests);

  Real computeValueLocal(const io_maps_type &x, const io_maps_type &params,
                         const io_maps_type &x_old) const override;
  io_maps_type
  computeVarGradientLocal(const io_maps_type &x, const io_maps_type &params,
                          const io_maps_type &x_old) const override;
  io_maps_type
  computeParamGradientLocal(const io_maps_type &x, const io_maps_type &params,
                            const io_maps_type &x_old) const override;

protected:
  Real dV1dqFun(const Real &q, const Real &Tn, const Real &D) const;
  virtual str_vct_type requiredParam() {
    str_vct_type temp = {"sVM", "sH", "eqedotc", "h", "n_V2"};
    return temp;
  };
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {"a"};
    return temp;
  };
};
