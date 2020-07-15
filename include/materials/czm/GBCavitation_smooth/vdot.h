#pragma once

#include "equationTests.h"
#include "nlFunRateNoHistory.h"
#include "pNormSoftMax.h"
#include "softMax.h"
#include "softMaxMin.h"
#include "v1dot.h"
#include "v2dotHigh.h"
#include "v2dotLow.h"

class vdot : public nlFunRateNoHistory {

public:
  vdot(const std::string &f_name, const v1dot &v1L, const v1dot &v1H,
       const v2dotLow &v2L, const v2dotHigh &v2H, const equationTests &eqTest,
       const int &max_type, const int &vdot_type,
       const int &triaxial_vdot_active, const Real &vdot_smooth_factor,
       const Real &theta_time_integration);

  // vdot(const std::string &fname, const Real &D, const Real &n, const Real
  // &h);

  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override;
  io_maps_type computeVarGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const override;
  io_maps_type computeParamGradient(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old) const override;

  const v1dot &_v1L;
  const v1dot &_v1H;
  const v2dotLow &_v2L;
  const v2dotHigh &_v2H;
  const softMaxMin _vmax;
  const equationTests &_eqTest;
  const pNormSoftMax _pnorm_max;
  const softMax _soft_max;
  const int _max_type;
  const int _vdot_type;
  const int _triaxial_vdot_active;
};
