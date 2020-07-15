#pragma once

#include "equationTests.h"
#include "nlFunRate.h"
#include "softMax.h"
#include "softMaxMin.h"

class bdot : public nlFunRate {

public:
  bdot(const std::string &f_name, const Real &beta, const Real &b0,
       const Real &bsat, const Real &Sthr, const Real &sigma0, const Real &FN,
       const equationTests &eqTest, const Real &theta_time_integration);

  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override;
  io_maps_type computeVarGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const override;
  io_maps_type computeParamGradient(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old) const override;

protected:
  const Real _beta;
  const Real _b0;
  const Real _bsat;
  const Real _Sthr;
  const Real _sigma0;
  const Real _FN;
  const equationTests &_eqTest;
  const softMaxMin _soft_max_min;
  const softMax _soft_max;
  const Real _TNmin = 1e-3;
};
