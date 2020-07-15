#pragma once

#include "nlFunRate.h"

class TdotBase : public nlFunRate {

public:
  TdotBase(const std::string &fname, const std::string &displacement_name,
           const std::string &traction_name, const Real &E,
           const Real &interface_thickness, const Real &theta_time_integration);

  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override;
  io_maps_type computeVarGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const override;
  io_maps_type computeParamGradient(const io_maps_type &x,
                                    const io_maps_type &params,
                                    const io_maps_type &x_old) const override;

protected:
  virtual Real Stiffness(const io_maps_type &x, const io_maps_type &params,
                         const io_maps_type &x_old) const;
  virtual io_maps_type StiffnessVarGradient(const io_maps_type &x,
                                            const io_maps_type &params,
                                            const io_maps_type &x_old) const;
  virtual TdotBase::io_maps_type
  StiffnessParamGradient(const io_maps_type &x, const io_maps_type &params,
                         const io_maps_type &x_old) const {
    UNUSED(x);
    UNUSED(params);
    UNUSED(x_old);
    return _empty_map;
  }

  const std::string _displacement_name;
  const std::string _displacement_dot_name;
  const std::string _traction_name;
  const Real _E;
  const Real _interface_thickness;
  virtual Real getSecondTermValue(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const = 0;
  virtual io_maps_type
  getSecondTermVarGradient(const io_maps_type &x, const io_maps_type &params,
                           const io_maps_type &x_old) const = 0;
  virtual io_maps_type
  getSecondTermParamGradient(const io_maps_type &x, const io_maps_type &params,
                             const io_maps_type &x_old) const = 0;
};
