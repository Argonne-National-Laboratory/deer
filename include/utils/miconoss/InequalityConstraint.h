#pragma once
#include "NLParameter.h"
#include "NLVar.h"
#include "miconosstype.h"

/**
Thie class implements inequality constraints. All constraints shall be
implemented as g<0
**/
class InequalityConstraint
{
public:
  InequalityConstraint(const uint lm_index,
                       const NLSystemVars & sys_vars,
                       const NLSystemParameters & sysparams,
                       const uint n_sys);
  /// method retuning the residual due to the
  vecD getR(const vecD & lm) const;
  vecD getJcolumn(const vecD & lm) const;
  vecD getJrow(const vecD & lm) const;
  bool constraintIsActive() const;
  uint getIndex() const { return _eq_index; };

protected:
  /// the actual constraint value
  virtual double gFun() const = 0;
  /// the derivative of the constraint, i.e. dg/dX
  virtual vecD dgFun_dx() const = 0;

  /// the langrance multiplier index this constraints refers to
  const uint _lm_index;
  /// the non linear system variables
  const NLSystemVars & _sys_vars;
  /// the non linear system parameters
  const NLSystemParameters & _sysparams;
  /// the number of variables excluding LM
  const uint _n_vars;
  /// the row in the nonlinear system
  const uint _eq_index;
  /// the total number of variables including LM
  const uint _n_sys;

private:
  /// method computing the actual value of the constraints
  double LM_equation(const vecD & lm) const;
};
