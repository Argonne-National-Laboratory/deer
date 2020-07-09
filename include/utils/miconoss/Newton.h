#pragma once

#include "NLParameter.h"
#include "NLSystem.h"
#include "miconosstype.h"

class Newton {
public:
  Newton(NLSystem *nlsys, NLSystemVars *_sys_vars,
         const double tolerance = 1e-6, const uint _max_iter = 50,
         const miconossmath::normtype normtype = miconossmath::normtype::L2);

  bool solve(vecD &lm, matrixD &J);
  bool solveSubstep(vecD &lm, matrixD &J, NLSystemParameters *const sysparams,
                    const std::vector<std::string> &pname, matrixD &Tangent,
                    bool &custom_interruption,
                    double &increment_at_custom_interruption,
                    const uint max_ncut = 1);

protected:
  void computeNewtonStep(const vecD &R, const matrixD &J);
  virtual bool customSubstepInterruption(NLSystemParameters *const sysparams) {
    return false;
  };
  NLSystem *const _nlsys;
  NLSystemVars *const _sys_vars;
  const double _tolerance;
  const uint _max_iter;
  const uint _n_eq;
  const uint _n_lm;
  const uint _n_total;
  vecD _dx;
  const miconossmath::normtype _normtype;
};
