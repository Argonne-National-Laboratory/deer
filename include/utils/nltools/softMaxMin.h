#pragma once

#include "nlFunBase.h"

class softMaxMin : public nlFunBase {

public:
  softMaxMin(const std::string & /*f_name*/, const Real &f = 50,
             const bool &max_of_abs = false,
             const bool &max_of_abs_signed = false,
             const bool &compute_min = false);

  // virtual void prepare(const io_maps_type &x,
  //                      const io_maps_type &params) override {
  //   updateX2useAndParam(x);
  // };

  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const;

  io_maps_type computeVarGradient(const io_maps_type &x,
                                  const io_maps_type &params,
                                  const io_maps_type &x_old) const override;
  io_maps_type
  computeParamGradient(const io_maps_type &x,
                       const io_maps_type &params,
                       const io_maps_type &x_old) const override {
    UNUSED(x);
    UNUSED(x_old);
    UNUSED(params);
    return _empty_map;
  };

  io_maps_type computeVarGradient(const io_maps_type & /*x*/,
                                  const map_of_io_maps_type & /*x_grad*/,
                                  const io_maps_type & /*params*/) const;
  io_maps_type computeParamGradient(const io_maps_type & /*x*/,
                                    const map_of_io_maps_type & /*x_grad*/,
                                    const io_maps_type & /*params*/) const;

protected:
  const Real _f;
  const bool _max_of_abs;
  const bool _max_of_abs_signed;
  const bool _compute_min;

  void updateParam(const io_maps_type & /*_x2use*/, Real & /*_c*/,
                   Real & /*_k*/, Real & /*_arg_log*/) const;
  void updateX2useAndParam(const io_maps_type & /*x*/,
                           io_maps_type & /*_x2use*/,
                           io_maps_type & /*_sign_x2use*/,
                           Real & /*_sign_xref*/) const;
  virtual str_vct_type requiredVar() {
    str_vct_type temp(0);
    return temp;
  };
  virtual str_vct_type requiredParam() {
    str_vct_type temp(0);
    return temp;
  };
};
