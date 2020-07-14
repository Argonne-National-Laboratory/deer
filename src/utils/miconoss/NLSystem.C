#include "NLSystem.h"

NLSystem::NLSystem(NLSystemVars *const sysvars,
                   NLSystemParameters *const sysparams,
                   std::vector<Equation *> eq,
                   std::vector<const InequalityConstraint *> lm,
                   NLPreEquationEvalautionCalc *pre_eval)
    : _sys_vars(sysvars), _sys_params(sysparams), _eq(eq), _lm(lm),
      _pre_eval(pre_eval), _nx(_eq.size()), _nlm(_lm.size()),
      _nsys(_nx + _nlm) {}

vecD NLSystem::assembleR(const vecD &lm) const {
  _pre_eval->updateValues();
  _pre_eval->updateDerivatives();

  vecD R(_nsys);

  for (uint i = 0; i < _nx; i++)
    R[_eq[i]->getIndex()] = _eq[i]->getR();

  for (uint j = 0; j < _nlm; j++) {
    vecD rlm = _lm[j]->getR(lm);
    for (uint i = 0; i < _nsys; i++)
      R[i] += rlm[i];
  }

  return R;
}

matrixD NLSystem::assembleJ(const vecD &lm) const {

  matrixD J(_nsys, vecD(_nsys, 0));

  for (uint i = 0; i < _nx; i++) {
    vecD jr = _eq[i]->getJrow();
    const uint eqidx = _eq[i]->getIndex();
    for (uint j = 0; j < _nx; j++)
      J[eqidx][j] = jr[j];
  }

  for (uint j = 0; j < _nlm; j++) {
    vecD jlm = _lm[j]->getJrow(lm);
    const uint jidx = _lm[j]->getIndex();

    for (uint k = 0; k < _nsys; k++)
      J[jidx][k] = jlm[k];

    jlm = _lm[j]->getJcolumn(lm);
    for (uint k = 0; k < _nx; k++)
      J[k][jidx] = jlm[k];
  }
  return J;
}

void NLSystem::updateEquationConstants() {
  _pre_eval->updateValuesExplicit();
  for (uint i = 0; i < _nx; i++)
    _eq[i]->updateConstants();
}

matrixD NLSystem::unscaleJacobian(const matrixD &Jscaled) {
  matrixD J = Jscaled;
  // scale equations
  for (uint i = 0; i < _nx; i++) {
    const double sf_x = _sys_vars->getDVarDVarScaled(i);
    for (uint j = 0; j < _nx; j++)
      J[i][j] = Jscaled[i][j] * sf_x * _sys_vars->getDVarScaledDVar(j);
  }

  // scale LM row nad columns equations
  for (uint i = _nx; i < _nsys; i++)
    for (uint j = 0; j < _nx; j++) {
      const double dx_scaled_dx = _sys_vars->getDVarScaledDVar(j);
      J[i][j] *= dx_scaled_dx;
      J[j][i] *= dx_scaled_dx;
    }

  return J;
}

matrixD NLSystem::getDResidualDParams(const std::vector<std::string> &pname) {
  updateEquationConstants();
  const uint n_param = pname.size();

  matrixD dresdp(n_param, vecD(_nsys));

  for (uint i = 0; i < _nx; i++) {
    _eq[i]->computeDEquationDP();
    for (uint p = 0; p < n_param; p++)
      dresdp[p][i] = _eq[i]->getDEquationDParam(pname[p]) *
                     _sys_vars->getDVarDVarScaled(i);
  }

  return dresdp;
}

matrixD NLSystem::getDSystemVarsDParams(const std::vector<std::string> &pname,
                                        const vecD &lm) {
  updateEquationConstants();
  matrixD J = assembleJ(lm);
  J = unscaleJacobian(J);
  const uint n_param = pname.size();

  matrixD dRdp = getDResidualDParams(pname);
  matrixD dxdP;
  int ierr = miconossmath::solveAxNb(J, dRdp, _nsys, dxdP);

  for (uint p = 0; p < n_param; p++)
    for (uint i = 0; i < _nx; i++)
      dxdP[p][i] *= -1;

  return dxdP;
}

void NLSystem::updateEquationScaling() const {
  for (uint i = 0; i < _nx; i++)
    _eq[i]->autoScaleEquation();
}
