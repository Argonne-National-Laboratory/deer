#include "InequalityConstraint.h"

InequalityConstraint::InequalityConstraint(const uint lm_index,
                                           const NLSystemVars &sys_vars,
                                           const NLSystemParameters &sysparams,
                                           const uint n_sys)
    : _lm_index(lm_index), _sys_vars(sys_vars), _sysparams(sysparams),
      _n_vars(_sys_vars.getNVars()), _eq_index(lm_index + _n_vars),
      _n_sys(n_sys) {}

vecD InequalityConstraint::getR(const vecD &lm) const {

  vecD R(_n_sys, 0);

  const vecD dgdx = dgFun_dx();
  for (uint i = 0; i < _n_vars; i++)
    R[i] = lm[_lm_index] * dgdx[i];

  R[_eq_index] = LM_equation(lm);
  return R;
}

vecD InequalityConstraint::getJcolumn(const vecD &lm) const {
  return dgFun_dx();
}

double InequalityConstraint::LM_equation(const vecD &lm) const {
  const double g_val = gFun();
  const double l_val = -lm[_lm_index];
  if (l_val < g_val)
    return g_val;
  else
    return l_val;
}

vecD InequalityConstraint::getJrow(const vecD &lm) const {
  const double g_val = gFun();
  const double l_val = -lm[_lm_index];
  vecD Jrow(_n_sys, 0);
  if (l_val < g_val) {
    const vecD dgdx = dgFun_dx();
    for (uint i = 0; i < _n_vars; i++)
      Jrow[i] = dgdx[i];
  } else
    Jrow[_eq_index] = -1;
  return Jrow;
}

bool InequalityConstraint::constraintIsActive() const { return gFun() > 0; }
