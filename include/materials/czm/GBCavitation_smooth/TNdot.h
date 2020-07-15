#pragma once

#include "TdotBase.h"
#include "vdot.h"

class TNdot : public TdotBase {

public:
  TNdot(const std::string &fname, const std::string &displacement_name,
        const std::string &traction_name, const Real &E,
        const Real &interface_thickness, const vdot &vDot,
        const Real &E_penalty, const Real &theta_time_integration);

protected:
  Real Stiffness(const io_maps_type &x, const io_maps_type &params,
                 const io_maps_type &x_old) const override;
  io_maps_type StiffnessVarGradient(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old) const override;

  bool innerPenetration(const Real &un) const;
  const vdot &_vdot;
  const Real _E_Penalty;

  Real getSecondTermValue(const io_maps_type &x, const io_maps_type &params,
                          const io_maps_type &x_old) const override;
  io_maps_type
  getSecondTermVarGradient(const io_maps_type &x, const io_maps_type &params,
                           const io_maps_type &x_old) const override;
  io_maps_type
  getSecondTermParamGradient(const io_maps_type &x, const io_maps_type &params,
                             const io_maps_type &x_old) const override;
};
