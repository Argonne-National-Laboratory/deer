#pragma once
#include "Equation.h"
#include "InequalityConstraint.h"
#include "NLPreEquationEvalautionCalc.h"
#include "miconosstype.h"

/** class representing a constraint or unconstrained non linear system.
It requires varaible (NLSystemVars), paramters (NLSystemParameters), a vector of
equations ( std::vector<Equation *>), a vector of lagrange multipliers (might be
empty),  and a pointer to a method computing intermedaite values used by
multiple equations (NLPreEquationEvalautionCalc). This class assumes teh user
properly setup indeces for varaible and equations, i.e. x_i =
eq_i(x1,...xN). THis calss is responsible for assembling the residual and the
jacobian
**/
class NLSystem {
public:
  NLSystem(NLSystemVars *const sysvars, NLSystemParameters *const sysparams,
           std::vector<Equation *> eq,
           const std::vector<const InequalityConstraint *> lm,
           NLPreEquationEvalautionCalc *pre_eval);
  vecD assembleR(const vecD &lm) const;
  matrixD assembleJ(const vecD &lm) const;

  /// method using NLPreEquationEvalautionCalc to precalcualte constants
  /// requried by the NLSystem
  void updateEquationConstants();
  /// method scaling equations
  void updateEquationScaling() const;

  /// get the number of NLvariables
  uint getNx() const { return _nx; };
  /// get the number of LM
  uint getNmultipliers() const { return _nlm; };
  /// get the NLSystem number of varaibles, inlcuding lagrange multipliers
  uint getDim() const { return _nlm + _nx; };

  /// get dResidual/dParam
  matrixD getDResidualDParams(const std::vector<std::string> &pname);

  /// get dNLVars/dParam
  matrixD getDSystemVarsDParams(const std::vector<std::string> &pname,
                                const vecD &lm);

  /// get the unscaled Jacobian
  matrixD unscaleJacobian(const matrixD &Jscaled);

protected:
  /// pointer to NLSystemVars
  NLSystemVars *const _sys_vars;
  /// pointer to NLSystemParameters
  NLSystemParameters *const _sys_params;
  /// vector of pointers to Equations
  std::vector<Equation *> _eq;
  /// vector of pointers to constraints
  std::vector<const InequalityConstraint *> _lm;
  /// pointer to precalculation method
  NLPreEquationEvalautionCalc *const _pre_eval;

  /// number of NLVariables
  const uint _nx;
  /// number of constraints
  const uint _nlm;
  /// total number of varaible plus constraints
  const uint _nsys;
};
