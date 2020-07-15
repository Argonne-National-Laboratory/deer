#pragma once

#include "bdot.h"
#include "nlFunConstant.h"
#include "nlFunRate.h"
#include "softMaxMin.h"

class bComp : public nlFunRate {

public:
  bComp(const std::string &fname, const bdot &bFun,
        const Real &theta_time_integration);

  // bComp(const std::string &f_name, const Real &h, const Real &n, const Real
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
  const bdot &_bdot;
  const softMaxMin _soft_min;
};
