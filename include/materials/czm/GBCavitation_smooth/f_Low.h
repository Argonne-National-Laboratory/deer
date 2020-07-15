#pragma once

#include "f_aL.h"
#include "f_ab.h"
#include "nlFunBase.h"
#include "softMax.h"
#include "softMaxMin.h"

class f_Low : public nlFunBase {

public:
  f_Low(const std::string &f_name, const f_aL &faL, const f_ab &fab);
  // f_Low(const std::string &f_name, const Real &D);
  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override;
  io_maps_type computeVarGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const override;
  io_maps_type computeParamGradient(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old) const override;

protected:
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {"a", "b"};
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp(0);
    return temp = {"eqedotc", "sVM"};
  };

  const f_aL &_faL;
  const f_ab &_fab;
  const softMaxMin _soft_max_min;
  const softMax _soft_max;
};
