#pragma once

#include "nlFunBase.h"

class qFun : public nlFunBase {

public:
  qFun(const std::string &f_name, const nlFunBase &f);

  // qFun(const std::string &fname);
  //
  // qFun(const std::string &fname, const Real &D);

  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override;
  io_maps_type computeVarGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const override;
  io_maps_type computeParamGradient(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old) const override;

protected:
  Real dqdf(const io_maps_type &x, const io_maps_type &params,
            const io_maps_type &x_old) const;
  virtual str_vct_type requiredVar() {
    str_vct_type temp(0);
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp(0);
    return temp;
  };

  const nlFunBase &_f;
};
