#pragma once
#include "NLParameter.h"
#include "NLVar.h"
#include "miconosstype.h"

class InequalityConstraint {
public:
  InequalityConstraint(const uint lm_index, const NLSystemVars &sys_vars,
                       const NLSystemParameters &sysparams, const uint n_sys);
  vecD getR(const vecD &lm) const;
  vecD getJcolumn(const vecD &lm) const;
  vecD getJrow(const vecD &lm) const;
  bool constraintIsActive() const;
  uint getIndex() const { return _eq_index; };

protected:
  /// implement gfun as constraint < 0
  virtual double gFun() const = 0;
  virtual vecD dgFun_dx() const = 0;

  const uint _lm_index;
  const NLSystemVars &_sys_vars;
  const NLSystemParameters &_sysparams;
  const uint _n_vars;
  const uint _eq_index;
  const uint _n_sys;

private:
  double LM_equation(const vecD &lm) const;
};
