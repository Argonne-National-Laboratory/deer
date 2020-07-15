//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef NLSOLVERVARTOOLS_H
#define NLSOLVERVARTOOLS_H

#include "Moose.h"

// Any requisite includes here
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

class var4NL {
public:
  var4NL(const std::vector<std::string> nls_var_names,
         const std::vector<std::string> other_var_names,
         const std::vector<std::string> parameter_names,
         const Real t_last_solution, const Real t_end);

  void resetAllVarValues();
  void resetAllVarDerivative();
  void resetVarDerivativeVector(const std::string & /*varname*/);
  unsigned int getVarIndex(const std::string & /*var_name*/) const;

  /// set time
  void setTimeInit(const Real &t_last_solution, const Real &t_end) {
    _t_last_solution = t_last_solution;
    _t_end = t_end;
    _t_current = t_end;
    _dt = _t_end - _t_last_solution;
  };

  void cutDt() {
    _dt = _dt / 2.;
    _t_current = _t_last_solution + _dt;
  };

  void reinitDt(const Real dt) {
    _dt = dt;
    _t_current = dt;
    _t_last_solution = 0;
  };

  void increaseDt() {
    _dt = _dt * 2.;
    _t_current = _t_last_solution + _dt;
    if (_t_current > _t_end) {
      _t_current = _t_end;
      _dt = _t_end - _t_last_solution;
    }
  };

  void advanceTime(bool &last_solution_computed) {
    last_solution_computed = false;

    if (_t_current == _t_end) {
      last_solution_computed = true;
      return;
    } else {
      _t_last_solution = _t_current;
      _t_current += _dt;
      if (_t_current > _t_end) {
        _t_current = _t_end;
        _dt = _t_current - _t_last_solution;
      }
    }
  };

  Real getDt() { return _dt; };
  Real getTcurrent() { return _t_current; };
  Real getTlastsolution() { return _t_last_solution; };

  /// SCALAR OPERATIONS:
  Real getVarValue(const std::string & /*v_name*/,
                   const bool scaled = false) const;
  void setVarValue(const std::string & /*v_name*/, const Real & /*v_value*/);

  void setVarMin(const std::string & /*v_name*/, const Real & /*v_value*/,
                 const DenseVector<Real> & /*dvmin_dx*/);

  void setVarMax(const std::string & /*v_name*/, const Real & /*v_value*/,
                 const DenseVector<Real> & /*dvmin_dx*/);
  Real getParameter(const std::string & /*p_name*/) const;
  void setParameter(const std::string & /*p_name*/, const Real & /*p_value*/);
  Real getVarLastSolution(const std::string & /*v_name*/) const;
  Real
  getVarLastSequentialSubstepSolution(const std::string & /*v_name*/) const;

  void setVarLastSolution(const std::string & /*v_name*/,
                          const Real & /*v_value*/);
  Real getVarResidual(const std::string & /*v_name*/) const;
  void setVarResidual(const std::string & /*v_name*/, const Real &
                      /*v_value*/);

  Real getDerivativeValue(const std::string & /*varname*/,
                          const std::string & /*dvarname*/) const;

  void setDerivativeValue(const std::string & /*varname*/,
                          const std::string & /*dvarname*/,
                          const Real & /*value*/);

  void setResidualDerivativeValue(const std::string & /*varname*/,
                                  const std::string & /*dvarname*/,
                                  const Real & /*value*/);

  DenseVector<Real> computeDx_DY(const std::string y_name) const;

  /// VECTOR OPERATIONS:
  DenseVector<Real>
  getDerivativeVectorValue(const std::string & /*varname*/) const;
  void setDerivativeVectorValue(const std::string & /*varname*/,
                                const DenseVector<Real> & /*d*/);

  void setResidualDerivativeVectorValue(const std::string & /*varname*/,
                                        const DenseVector<Real> & /*d*/);
  /// copy stuff
  DenseVector<Real> getAllValue() const { return _var_value; };
  DenseVector<Real> getAllResidual() const { return _var_residual; };
  DenseMatrix<Real> getAllDerivatives() const { return _var_derivatives; };
  DenseVector<Real> getAllParameters() const { return _parameters; };

  void getALL(DenseVector<Real> &v, DenseVector<Real> &r, DenseMatrix<Real> &D,
              DenseVector<Real> &p) const {
    v = getAllValue();
    r = getAllResidual();
    D = getAllDerivatives();
    p = getAllParameters();
  }

  void setAllValue(const DenseVector<Real> &values) {
    for (unsigned int i = 0; i < _n_all_var_names; i++)
      _var_value(i) = values(i);
  };

  void setAllResidual(const DenseVector<Real> &values) {
    for (unsigned int i = 0; i < _n_all_var_names; i++)
      _var_residual(i) = values(i);
  };

  void setAllDerivatives(const DenseMatrix<Real> &values) {
    for (unsigned int i = 0; i < _n_all_var_names; i++)
      for (unsigned int j = 0; j < _n_all_var_names; j++)
        _var_derivatives(i, j) = values(i, j);
  };

  void setAllParameters(const DenseVector<Real> &values) {
    for (unsigned int i = 0; i < _n_all_var_names; i++)
      _parameters(i) = values(i);
  };

  void setALL(const DenseVector<Real> &v, const DenseVector<Real> &r,
              const DenseMatrix<Real> &D, const DenseVector<Real> &p) {
    setAllValue(v);
    setAllResidual(r);
    setAllDerivatives(D);
    setAllParameters(p);
  }

  const Real _dvar_fdiff = 1e-6;
  Real perturb(const std::string &name) {
    unsigned int idx = getVarIndex(name);
    Real v_pert = 0;
    Real delta = 0;
    if (idx < (_n_all_var_names - _n_parameters))
      v_pert = getVarValue(name);
    else
      v_pert = getParameter(name);

    if (std::abs(v_pert) > 0.)
      delta = v_pert * _dvar_fdiff;
    else
      delta = _dvar_fdiff;

    if (std::copysign(1., v_pert) != std::copysign(1., v_pert + delta))
      mooseError("perturbation changed the variable sign");
    v_pert += delta;

    if (idx < (_n_all_var_names - _n_parameters))
      setVarValue(name, v_pert);
    else
      setParameter(name, v_pert);

    return delta;
  }
  // #endif
  /// floating vectors
  DenseVector<Real> getNewZeroDerivativeVector() const;
  DenseVector<Real>
  getNewUnitDerivativeVector(const std::string & /*v_name*/) const;
  void setZeroDerivativeVectorValue(DenseVector<Real> & /*dvector*/,
                                    const std::string & /*dvarname*/,
                                    const Real & /*value*/);

  /// NON LINEAR SYSTEM OPERATIONS
  DenseVector<Real> getValueNL(const bool scaled = false) const;
  void setValueNL(const DenseVector<Real> & /*value*/,
                  const bool from_scaled = false);
  DenseVector<Real> getResidualNL(const bool scaled = false) const;
  void setResidualNL(DenseVector<Real> &residual);
  DenseVector<Real> getLastSolutionNL() const;
  DenseMatrix<Real> getJacobianNL(const bool scaled = false) const;
  void setScaleFactorNL();

  DenseVector<Real> getScaleFactorNL();
  void checkALL() const;
  void resetToLastSolution();
  void advanceStep();
  void printNLTimes() const;
  void printNLVars() const;
  void printNLDerivatives() const;
  void setVarScaleFactor(const std::string &vname, const Real &sf) {
    unsigned int v_idx = getVarIndex(vname);
    _scale_factor[v_idx] = sf;
    if (_scale_factor[v_idx] == 0.)
      _scale_factor[v_idx] = 1.;
  }
  Real getVarScaleFactor(const std::string &vname) {
    unsigned int v_idx = getVarIndex(vname);
    return _scale_factor[v_idx];
  }
  bool check_vars_are_finite() const;
  bool check_derivatives_are_finite() const;
  void SequentialSubstepOn() { _sequential_substep = true; };
  void SequentialSubstepOff() { _sequential_substep = false; };

protected:
  const unsigned int _n_nls_vars;
  const std::vector<std::string> _nls_var_names;
  const unsigned int _n_other_vars;
  const std::vector<std::string> _other_var_names;
  const unsigned int _n_parameters;
  const std::vector<std::string> _parameter_names;
  const unsigned int _n_all_var_names;
  const std::vector<std::string> _all_var_names;
  const std::map<std::string, unsigned int> _map_var_name_idx;
  DenseVector<Real> _var_value;
  DenseVector<Real> _var_residual;
  DenseMatrix<Real> _var_derivatives;
  DenseMatrix<Real> _var_residual_derivatives;
  DenseVector<Real> _var_last_solution;
  DenseVector<Real> _var_last_sequential_substep_solution;
  DenseVector<Real> _parameters;
  Real _t_last_solution;
  Real _t_end;
  Real _t_current;
  Real _dt;
  std::vector<Real> _scale_factor;
  Real _alpha;

  /// concatenate vector string
  std::vector<std::string>
  InitAllVarVarNames(const std::vector<std::vector<std::string>> /*v*/) const;

  /// initialize a mab between var names and indices
  std::map<std::string, unsigned int>
  InitVarIdxMap(const std::vector<std::string> /* var_names */) const;
  bool _sequential_substep = false;
};

#endif // NLSOLVERVARTOOLS_H
