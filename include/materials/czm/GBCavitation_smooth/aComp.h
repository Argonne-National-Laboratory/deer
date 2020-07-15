#pragma once

#include "adot.h"
#include "nlFunConstant.h"
#include "nlFunRate.h"
#include "softMax.h"
#include "softMaxMin.h"

class aComp : public nlFunRate {

public:
  aComp(const std::string &fname, const adot &aFun, const Real &amin,
        const Real &theta_time_integration);

  // aComp(const std::string &f_name, const Real &h, const Real &n, const Real
  // &D);

  Real computeTimeIntegral(const io_maps_type &x, const io_maps_type &params,
                           const io_maps_type &x_old,
                           const Real &dt) const override;

  io_maps_type computeTimeIntegralVarGradient(const io_maps_type &x,
                                              const io_maps_type &params,
                                              const io_maps_type &x_old,
                                              const Real &dt) const override;

  io_maps_type computeTimeIntegralParamGradient(const io_maps_type &x,
                                                const io_maps_type &params,
                                                const io_maps_type &x_old,
                                                const Real &dt) const override;

protected:
  const adot &_adot;
  const nlFunConstant _a0Fun;
  const softMaxMin _soft_max_min;
  const softMax _soft_max;
};
