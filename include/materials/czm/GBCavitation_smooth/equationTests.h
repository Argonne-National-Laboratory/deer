#pragma once

#include "nlFunBase.h"

class equationTests : public nlFunBase {

public:
  equationTests(const std::string &f_name, const Real &a0, const Real &b0,
                const Real &sigma_0, const Real &S_thr,
                const Real &beta_exponent, const Real &b_saturation,
                const Real &E_interface, const bool &cavity_nucleation_on,
                const bool &cavity_growth_on,
                const bool &triaxial_cavity_growth_on);

  bool nucleationIsActive(const io_maps_type &x, const io_maps_type &params,
                          const io_maps_type &x_old) const;

  bool growthIsActive(const io_maps_type &x, const io_maps_type &params,
                      const io_maps_type &x_old) const;

  bool triaxialGrowthIsActive(const io_maps_type &x, const io_maps_type &params,
                              const io_maps_type &x_old) const;

  bool linearUpdate(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const;

  bool nucleationAboveThreshold(const io_maps_type &params,
                                const io_maps_type &x_old) const;

protected:
  const Real _a0;
  const Real _b0;
  const Real _sigma_0;
  const Real _S_thr;
  const Real _beta_exponent;
  const Real _b_saturation;
  const Real _E_interface;
  const bool _cavity_nucleation_on;
  const bool _cavity_growth_on;
  const bool _triaxial_cavity_growth_on;

  virtual str_vct_type requiredVar() {
    str_vct_type temp = {"a", "b"};
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp(0);
    return temp;
  };
};
