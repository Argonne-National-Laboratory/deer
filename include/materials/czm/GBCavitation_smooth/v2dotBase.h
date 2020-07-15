#pragma once

#include "equationTests.h"
#include "logistic.h"
#include "nlFunRateNoHistory.h"
#include "softMax.h"

class v2dotBase : public nlFunRateNoHistory {

public:
  v2dotBase(const std::string &f_name, const Real &n, const Real &h,
            const Real &theta_time_integration, const equationTests &eqTests);

  Real computeValue(const io_maps_type &x, const io_maps_type &params,
                    const io_maps_type &x_old) const override final {

    if (!_eqTests.triaxialGrowthIsActive(x, params, x_old))
      return 0;

    return computeValueLocal(x, params, x_old);

    // Real TN = getValueFromMap(x, "Tn", "v2dotBase::computeValue");
    // Real logitTN = _logit.computeLogisticValue(TN);
    // return computeValueLocal(x, params, x_old) * logitTN;
  };

  io_maps_type
  computeVarGradient(const io_maps_type &x, const io_maps_type &params,
                     const io_maps_type &x_old) const override final {

    if (!_eqTests.triaxialGrowthIsActive(x, params, x_old))
      return _empty_map;

    return computeVarGradientLocal(x, params, x_old);

    // Real TN = getValueFromMap(x, "Tn", "v2dotBase::computeValue");
    // Real g = _logit.computeLogisticValue(TN);
    // io_maps_type dlogit_dtn = _logit.computeLogisticGradient(TN);
    // io_maps_type dg_dx = {{"Tn", dlogit_dtn.find("x")->second}};
    //
    // Real f = computeValueLocal(x, params, x_old);
    // io_maps_type df_dx = computeVarGradientLocal(x, params, x_old);
    //
    // return f_times_g_D(f, df_dx, g, dg_dx);
  };

  io_maps_type
  computeParamGradient(const io_maps_type &x, const io_maps_type &params,
                       const io_maps_type &x_old) const override final {
    if (!_eqTests.triaxialGrowthIsActive(x, params, x_old))
      return _empty_map;
    return computeParamGradientLocal(x, params, x_old);
  };

protected:
  void commonOperation(const io_maps_type &x, const io_maps_type &params,
                       Real &sVM, Real &sH, Real &triax, Real &e, Real &a,
                       Real &alpha, Real &m, Real &beta) const;

  Real dV1dqFun(const Real &q, const Real &Tn, const Real &D) const;
  virtual str_vct_type requiredParam() {
    str_vct_type temp = {"sVM", "sH", "eqedotc", "h", "n_V2"};
    return temp;
  };
  virtual str_vct_type requiredVar() {
    str_vct_type temp = {"a"};
    return temp;
  };
  const Real _n;
  const Real _h;
  const Real _alpha;

  virtual Real computeValueLocal(const io_maps_type &x,
                                 const io_maps_type &params,
                                 const io_maps_type &x_old) const = 0;

  virtual io_maps_type
  computeVarGradientLocal(const io_maps_type &x, const io_maps_type &params,
                          const io_maps_type &x_old) const = 0;
  virtual io_maps_type
  computeParamGradientLocal(const io_maps_type &x, const io_maps_type &params,
                            const io_maps_type &x_old) const = 0;

  Real alphaFun(const Real &n) const;
  Real gFun(const Real &m) const;
  Real betaFun(const Real &m) const;
  const equationTests &_eqTests;
  const logistic _logit;
};
