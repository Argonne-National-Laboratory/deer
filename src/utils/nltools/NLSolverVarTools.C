//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseError.h"
#include "NLSolverVarTools.h"

var4NL::var4NL(const std::vector<std::string> nls_var_names,
               const std::vector<std::string> other_var_names,
               const std::vector<std::string> parameter_names,
               const Real t_last_solution, const Real t_end)
    : _n_nls_vars(nls_var_names.size()), _nls_var_names(nls_var_names),
      _n_other_vars(other_var_names.size()), _other_var_names(other_var_names),
      _n_parameters(parameter_names.size()), _parameter_names(parameter_names),
      _n_all_var_names(_n_nls_vars + _n_other_vars + _n_parameters),
      _all_var_names(InitAllVarVarNames(
          {_nls_var_names, other_var_names, parameter_names})),
      _map_var_name_idx(InitVarIdxMap(_all_var_names)),
      _var_value(_n_all_var_names), _var_residual(_n_all_var_names),
      _var_derivatives(_n_all_var_names, _n_all_var_names),
      _var_residual_derivatives(_n_all_var_names, _n_all_var_names),
      _var_last_solution(_n_all_var_names),
      _var_last_sequential_substep_solution(_n_all_var_names),
      _parameters(_n_all_var_names), _t_last_solution(t_last_solution),
      _t_end(t_end), _t_current(t_last_solution),
      _dt(_t_last_solution - _t_end), _scale_factor(_n_nls_vars, 1.)

{}

std::vector<std::string> var4NL::InitAllVarVarNames(
    const std::vector<std::vector<std::string>> v) const {
  unsigned int i, j, n, m;
  std::vector<std::string> t(0);
  n = v.size();

  for (i = 0; i < n; i++) {
    m = v[i].size();
    for (j = 0; j < m; j++)
      t.push_back(v[i][j]);
  }
  return t;
}

std::map<std::string, unsigned int>
var4NL::InitVarIdxMap(const std::vector<std::string> var_names) const {
  std::map<std::string, unsigned int> map_var_name_idx;
  map_var_name_idx.clear();
  for (unsigned int i = 0; i < var_names.size(); i++) {
    map_var_name_idx[var_names[i]] = i;
  }
  return map_var_name_idx;
}

void var4NL::resetVarDerivativeVector(const std::string &varname) {
  unsigned int v_idx = getVarIndex(varname);
  for (unsigned int i = 0; i < _n_all_var_names; i++) {
    _var_derivatives(v_idx, i) = 0;
    if (i == v_idx)
      _var_derivatives(v_idx, i) = +1;
  }
}

unsigned int var4NL::getVarIndex(const std::string &var_name) const {
  auto var_idx_p = _map_var_name_idx.find(var_name);
  if (var_idx_p != _map_var_name_idx.end())
    return var_idx_p->second;
  else
    mooseError("getVarIndex cannot find variable named ", var_name);
}

void var4NL::resetAllVarDerivative() {
  for (unsigned int i = 0; i < _n_all_var_names; i++)
    for (unsigned int j = 0; j < _n_all_var_names; j++) {
      _var_derivatives(i, j) = 0;
      if (i == j)
        _var_derivatives(i, j) += 1;
    }
}

void var4NL::resetAllVarValues() {
  for (unsigned int i = 0; i < _n_all_var_names; i++) {
    _var_value(i) = 0;
    _var_residual(i) = 0;
  }
}

void var4NL::resetToLastSolution() {
  resetAllVarDerivative();
  resetAllVarValues();
  for (unsigned int i = 0; i < _n_all_var_names; i++)
    if (!_sequential_substep)
      _var_value(i) = _var_last_solution(i);
    else
      _var_value(i) = _var_last_sequential_substep_solution(i);
}

void var4NL::advanceStep() {
  if (!_sequential_substep) {
    for (unsigned int i = 0; i < _n_all_var_names; i++) {
      _var_last_solution(i) = 0;
      if (i < _n_nls_vars)
        _var_last_solution(i) = _var_value(i);
    }
  } else {
    for (unsigned int i = 0; i < _n_all_var_names; i++) {
      if (i < _n_nls_vars)
        _var_last_sequential_substep_solution(i) = _var_value(i);
    }
    resetToLastSolution();
  }
}

void var4NL::setZeroDerivativeVectorValue(DenseVector<Real> &dvector,
                                          const std::string &dvarname,
                                          const Real &value) {
  dvector(getVarIndex(dvarname)) = value;
}

DenseVector<Real> var4NL::getNewZeroDerivativeVector() const {
  DenseVector<Real> d_vector(_n_all_var_names);
  return d_vector;
}

DenseVector<Real>
var4NL::getNewUnitDerivativeVector(const std::string &v_name) const {
  unsigned int v_idx = getVarIndex(v_name);
  DenseVector<Real> d_vector = getNewZeroDerivativeVector();
  d_vector(v_idx) = 1;
  return d_vector;
}

Real var4NL::getVarValue(const std::string &v_name, const bool scaled) const {
  const unsigned int v_idx = getVarIndex(v_name);
  if (!scaled)
    return _var_value(v_idx);
  else
    return _var_value(getVarIndex(v_name)) / _scale_factor[getVarIndex(v_name)];
}

void var4NL::setVarValue(const std::string &v_name, const Real &v_value) {
  _var_value(getVarIndex(v_name)) = v_value;
}

Real var4NL::getParameter(const std::string &p_name) const {
  return _parameters(getVarIndex(p_name));
}
void var4NL::setParameter(const std::string &p_name, const Real &p_value) {
  _parameters(getVarIndex(p_name)) = p_value;
}

Real var4NL::getVarLastSolution(const std::string &v_name) const {
  return _var_last_solution(getVarIndex(v_name));
}

Real var4NL::getVarLastSequentialSubstepSolution(
    const std::string &v_name) const {
  return _var_last_sequential_substep_solution(getVarIndex(v_name));
}

void var4NL::setVarLastSolution(const std::string &v_name,
                                const Real &v_value) {
  _var_last_solution(getVarIndex(v_name)) = v_value;
  if (_sequential_substep)
    _var_last_sequential_substep_solution(getVarIndex(v_name)) = v_value;
}

DenseVector<Real> var4NL::getLastSolutionNL() const {
  DenseVector<Real> d(_n_nls_vars);
  unsigned int i;
  for (i = 0; i < _nls_var_names.size(); i++)
    d(i) = getVarLastSolution(_nls_var_names[i]);
  return d;
}

DenseVector<Real> var4NL::getValueNL(const bool scaled) const {
  DenseVector<Real> d(_n_nls_vars);
  unsigned int i;
  for (i = 0; i < _n_nls_vars; i++)
    d(i) = getVarValue(_nls_var_names[i], scaled);
  return d;
}

void var4NL::setValueNL(const DenseVector<Real> &value,
                        const bool from_scaled) {
  DenseVector<Real> d(_n_nls_vars);
  unsigned int i;

  for (i = 0; i < _n_nls_vars; i++)
    if (from_scaled)
      setVarValue(_nls_var_names[i], value(i) * _scale_factor[i]);
    else
      setVarValue(_nls_var_names[i], value(i));
}

DenseMatrix<Real> var4NL::getJacobianNL(const bool scaled) const {
  DenseMatrix<Real> jac(_n_nls_vars, _n_nls_vars);
  unsigned int i, j;
  for (i = 0; i < _n_nls_vars; i++)
    for (j = 0; j < _n_nls_vars; j++) {
      jac(i, j) = _var_residual_derivatives(i, j);
      if (scaled)
        jac(i, j) *= _scale_factor[j] / _scale_factor[i];
    }
  return jac;
}

void var4NL::setScaleFactorNL() {
  unsigned int i;
  for (i = 0; i < _n_nls_vars; i++) {
    Real sf = std::abs(_var_last_solution(i));
    if (sf == 0.)
      sf = 1.;
    _scale_factor[i] = sf;
  }
}

DenseVector<Real> var4NL::getScaleFactorNL() {
  DenseVector<Real> sf(_n_nls_vars);
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    sf(i) = _scale_factor[i];
  return sf;
}

Real var4NL::getVarResidual(const std::string &v_name) const {
  unsigned int v_idx = getVarIndex(v_name);
  if (v_idx < _n_nls_vars)
    return _var_residual(v_idx);
  else
    mooseError("::getVarResidual only non linear variables has "
               "a Residual");
}

DenseVector<Real> var4NL::getResidualNL(const bool scaled) const {
  DenseVector<Real> d(_n_nls_vars);
  unsigned int i;
  for (i = 0; i < _n_nls_vars; i++) {
    d(i) = getVarResidual(_nls_var_names[i]);
    if (scaled)
      d(i) /= _scale_factor[i];
  }
  return d;
}

void var4NL::setResidualNL(DenseVector<Real> &residual) {
  unsigned int i;
  for (i = 0; i < _n_nls_vars; i++)
    setVarResidual(_nls_var_names[i], residual(i));
}

void var4NL::setVarResidual(const std::string &v_name, const Real &v_value) {
  unsigned int v_idx = getVarIndex(v_name);
  if (v_idx < _n_nls_vars)
    _var_residual(v_idx) = v_value;
  else
    mooseError("::setVarResidual only non linear variables has "
               "a Residual");
}

Real var4NL::getDerivativeValue(const std::string &varname,
                                const std::string &dvarname) const {
  return _var_derivatives(getVarIndex(varname), getVarIndex(dvarname));
}

DenseVector<Real>
var4NL::getDerivativeVectorValue(const std::string &varname) const {
  DenseVector<Real> d(_n_all_var_names);
  unsigned int v_idx = getVarIndex(varname);

  for (unsigned int i = 0; i < _n_all_var_names; i++)
    d(i) = _var_derivatives(v_idx, i);
  return d;
}

void var4NL::setDerivativeValue(const std::string &varname,
                                const std::string &dvarname,
                                const Real &value) {
  _var_derivatives(getVarIndex(varname), getVarIndex(dvarname)) = value;
}

void var4NL::setResidualDerivativeValue(const std::string &varname,
                                        const std::string &dvarname,
                                        const Real &value) {
  _var_residual_derivatives(getVarIndex(varname), getVarIndex(dvarname)) =
      value;
}

void var4NL::setDerivativeVectorValue(const std::string &varname,
                                      const DenseVector<Real> &d) {
  if (d.size() == _n_all_var_names) {
    unsigned int v_idx = getVarIndex(varname);
    for (unsigned int i = 0; i < _n_all_var_names; i++)
      _var_derivatives(v_idx, i) = d(i);
  } else
    mooseError("setDerivativeVectorValue: the size of the input "
               "vector does not match the number of variables");
}

void var4NL::setResidualDerivativeVectorValue(const std::string &varname,
                                              const DenseVector<Real> &d) {
  if (d.size() == _n_all_var_names) {
    unsigned int v_idx = getVarIndex(varname);
    for (unsigned int i = 0; i < _n_all_var_names; i++)
      _var_residual_derivatives(v_idx, i) = d(i);
  } else
    mooseError("setResidualDerivativeVectorValue: the size of the input "
               "vector does not match the number of variables");
}

bool var4NL::check_vars_are_finite() const {
  bool check = true;
  for (unsigned int i = 0; i < _n_all_var_names; i++)
    if (!std::isfinite(_var_value(i)) || !std::isfinite(_var_residual(i))) {
      check = false;
      break;
    }
  return check;
}

bool var4NL::check_derivatives_are_finite() const {
  bool check = true;
  for (unsigned int i = 0; i < _n_all_var_names; i++) {
    for (unsigned int j = 0; j < _n_all_var_names; j++)
      if (!std::isfinite(_var_derivatives(i, j)) ||
          !std::isfinite(_var_residual_derivatives(i, j))) {
        check = false;
        break;
      }
    if (!check)
      break;
  }
  return check;
}

void var4NL::printNLTimes() const {
  unsigned int i;

  std::cout << " NL TIMES   " << std::endl;

  std::cout << "_t_last_solution: " << _t_last_solution << std::endl;
  std::cout << "_t_current: " << _t_current << std::endl;
  std::cout << "_t_end: " << _t_end << std::endl;
  std::cout << "_dt: " << _dt << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

void var4NL::printNLVars() const {
  unsigned int i;

  std::cout << "   VARIABLE VALUES   " << std::endl;
  for (i = 0; i < _n_all_var_names - _n_parameters; i++) {
    std::cout << _all_var_names[i] << ": " << _var_value(i);
    if (i < _n_nls_vars) {
      std::cout << "   last solutionvalue: " << _var_last_solution(i);
      std::cout << "   scale factor: " << _scale_factor[i];
    }
    std::cout << std::endl;
  }
  std::cout << "   PARAMETERS VALUES   " << std::endl;
  for (i = _n_all_var_names - _n_parameters; i < _n_all_var_names; i++) {
    std::cout << _all_var_names[i] << ": " << _parameters(i) << std::endl;
  }
  std::cout << "   RESIDUAL VALUES   " << std::endl;
  for (i = 0; i < _n_nls_vars; i++) {
    std::cout << _all_var_names[i] << ": " << _var_residual(i);
    std::cout << std::endl;
  }
}

void var4NL::printNLDerivatives() const {
  unsigned int i, j;

  std::cout << "   VARIABLE DERIVATIVES   " << std::endl;
  for (i = 0; i < _n_all_var_names; i++)
    std::cout << " d" << _all_var_names[i] << " ";
  std::cout << std::endl;
  for (i = 0; i < _n_all_var_names; i++) {
    for (j = 0; j < _n_all_var_names; j++) {
      if (j == 0)
        std::cout << "d" << _all_var_names[i] << " ";
      std::cout << _var_derivatives(i, j) << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "   RESIDUAL DERIVATIVES   " << std::endl;
  for (i = 0; i < _n_all_var_names; i++)
    std::cout << " d" << _all_var_names[i] << " ";
  std::cout << std::endl;
  for (i = 0; i < _n_nls_vars; i++) {
    for (j = 0; j < _n_all_var_names; j++) {
      if (j == 0)
        std::cout << "d" << _all_var_names[i] << " ";
      std::cout << _var_residual_derivatives(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

void var4NL::checkALL() const {
  bool check_var = check_vars_are_finite();
  bool check_der = check_derivatives_are_finite();
  bool check_all = check_var && check_der;
  if (!check_all) {
    printNLVars();
    printNLDerivatives();
    mooseError("some of the variables or derivatives are not finite");
  }
}

DenseVector<Real> var4NL::computeDx_DY(const std::string y_name) const {
  DenseVector<Real> dR_dy(_n_nls_vars);
  DenseVector<Real> dx_dy(_n_nls_vars);
  DenseMatrix<Real> J = getJacobianNL(/*from_scaled =*/false);
  for (unsigned int i = 0; i < _n_nls_vars; i++) {
    std::string v_name = _nls_var_names[i];
    dR_dy(i) = getDerivativeValue(v_name + "_comp", y_name);
  }
  J.lu_solve(dR_dy, dx_dy);
  return dx_dy;
}
