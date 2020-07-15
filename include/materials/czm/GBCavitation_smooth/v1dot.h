#pragma once

#include "nlFunRateNoHistory.h"
#include "qFun.h"

class v1dot : public nlFunRateNoHistory {

public:
  v1dot(const std::string &f_name, const qFun &q, const Real &D,
        const Real &theta_time_integration);
  // v1dot(const std::string &fname, const Real &D, const bool &high,
  //       const Real &theta_time_integration);

  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override;
  io_maps_type computeVarGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const override;
  io_maps_type computeParamGradient(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old) const override;

protected:
  Real dV1dqFun(const Real &q, const Real &Tn) const;
  virtual str_vct_type requiredParam() {
    str_vct_type temp;
    return temp;
  };
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {"Tn"};
    return temp;
  };
  const qFun &_q;
  const Real _D;
};
