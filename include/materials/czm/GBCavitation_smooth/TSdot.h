#pragma once

#include "TdotBase.h"

class TSdot : public TdotBase {

public:
  TSdot(const std::string &fname, const std::string &displacement_name,
        const std::string &traction_name, const Real &E,
        const Real &interface_thickness, const Real &eta,
        const Real &theta_time_integration);

protected:
  const Real _eta;

  Real getSecondTermValue(const io_maps_type &x, const io_maps_type &params,
                          const io_maps_type &x_old) const override;
  io_maps_type
  getSecondTermVarGradient(const io_maps_type &x, const io_maps_type &params,
                           const io_maps_type &x_old) const override;
  io_maps_type
  getSecondTermParamGradient(const io_maps_type &x, const io_maps_type &params,
                             const io_maps_type &x_old) const override;
};
