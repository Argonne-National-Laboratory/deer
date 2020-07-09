#include "Equation.h"

#include <iostream>

Equation::Equation(const uint eq_index, NLSystemVars &sysvars,
                   NLSystemParameters &sysparams,
                   NLPreEquationEvalautionCalc &pre_eval)
    : _eq_index(eq_index), _sys_vars(sysvars), _sysparams(sysparams),
      _n_vars(_sys_vars.getNVars()), _n_params(_sysparams.getNParams()),
      _pre_eval(pre_eval), _dequation_dparam(_n_params) {}

void Equation::checkGradinet(double tol, double eps) {
  _pre_eval.updateAll(true);
  updateConstants();
  const double R = getR();
  const vecD dR_dx = getJrow();

  for (uint i = 0; i < _n_vars; i++) {
    const double old_var_value = _sys_vars.getValue(i);
    _sys_vars.setValue(i, old_var_value + eps);
    _pre_eval.updateAll(true);
    updateConstants();
    const double Rcurrent = getR();
    const double deltaR = Rcurrent - R;
    const double dR_dxi = (deltaR) / eps;
    const double diff = std::abs(dR_dxi - dR_dx[i]);
    if (diff > tol) {
      std::cout << "Equation " << _eq_index
                << ": FD and analytial derivatives do not match: \n FD = "
                << dR_dxi << " an= " << dR_dx[i] << " diff= " << diff
                << " component " << i << ". Delta R: " << deltaR << "\n\n";
      std::cout << "current R " << Rcurrent << " initial R " << R << "\n";
    }
    _sys_vars.setValue(i, old_var_value);
  }
}

RateEquation::RateEquation(const uint eq_index, NLSystemVars &sysvars,
                           NLSystemParameters &sysparams,
                           NLPreEquationEvalautionCalc &pre_eval,
                           const double theta)
    : Equation(eq_index, sysvars, sysparams, pre_eval), _theta(theta),
      _dexplicit_rate_dp(_n_params) {}

void RateEquation::updateConstants() {
  if (_theta != 0) {
    _explicit_rate = computedRate(/*implicit = */ false);
    _dexplicit_rate_dp = DComputedRatetDP(/*implicit =*/false);
  }
}

double RateEquation::getR() const {
  return _sys_vars.getValueScaled(_eq_index) -
         computedVal() / _sys_vars.getScalingFactor(_eq_index);
}

vecD RateEquation::getJrow() const {
  vecD Rgrad(_n_vars, 0);
  Rgrad[_eq_index] = 1;
  const double sf_x = _sys_vars.getScalingFactor(_eq_index);
  vecD dCompVal_dxi = DComputedValDx();
  for (uint i = 0; i < _n_vars; i++)
    Rgrad[i] -= dCompVal_dxi[i] / sf_x * _sys_vars.getDVarDVarScaled(i);
  return Rgrad;
}

double RateEquation::computedVal() const {

  return _sys_vars.getValueOld(_eq_index) +
         (_theta * _explicit_rate +
          (1. - _theta) * computedRate(/*implicit =*/true)) *
             _sysparams.getValue("dt");
}

vecD RateEquation::DComputedValDx() const {
  vecD dInc_Dx = DComputedRatetDx(/*implicit =*/true);
  for (uint i = 0; i < _n_vars; i++)
    dInc_Dx[i] *= (1. - _theta) * _sysparams.getValue("dt");
  return dInc_Dx;
}

void RateEquation::computeDEquationDP() {
  vecD dInc_Dp = DComputedRatetDP(/*implicit =*/true);
  _dexplicit_rate_dp = DComputedRatetDP(/*implicit =*/false);

  const double sf_x = _sys_vars.getDVarScaledDVar(_eq_index);
  for (uint i = 0; i < _n_params; i++)
    _dequation_dparam[i] =
        -(_theta * _dexplicit_rate_dp[i] + (1. - _theta) * dInc_Dp[i]) *
        _sysparams.getValue("dt") * sf_x;
}
