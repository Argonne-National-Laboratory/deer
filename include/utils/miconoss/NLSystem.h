#pragma once
#include "Equation.h"
#include "InequalityConstraint.h"
#include "NLPreEquationEvalautionCalc.h"
#include "miconosstype.h"

class NLSystem {
public:
  NLSystem(NLSystemVars *const sysvars, NLSystemParameters *const sysparams,
           std::vector<Equation *> eq,
           const std::vector<const InequalityConstraint *> lm,
           NLPreEquationEvalautionCalc *pre_eval);
  vecD assembleR(const vecD &lm) const;
  matrixD assembleJ(const vecD &lm) const;

  void updateEquationConstants();
  void updateEquationScaling() const;

  uint getNx() const { return _nx; };
  uint getNmultipliers() const { return _nlm; };
  uint getDim() const { return _nlm + _nx; };

  matrixD getDResidualDParams(const std::vector<std::string> &pname);

  matrixD getDSystemVarsDParams(const std::vector<std::string> &pname,
                                const vecD &lm);

  matrixD unscaleJacobian(const matrixD &Jscaled);

protected:
  NLSystemVars *const _sys_vars;
  NLSystemParameters *const _sys_params;
  std::vector<Equation *> _eq;
  std::vector<const InequalityConstraint *> _lm;
  NLPreEquationEvalautionCalc *const _pre_eval;

  const uint _nx;
  const uint _nlm;
  const uint _nsys;
};
