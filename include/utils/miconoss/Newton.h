#pragma once

#include "NLParameter.h"
#include "NLSystem.h"
#include "miconosstype.h"

/*
this class solves a constrained or unconstrained nonlinear sytem using newton
method. To solve the system the user shall call solve or solveSubstep
*/
class Newton
{
public:
  Newton(NLSystem * nlsys,
         NLSystemVars * _sys_vars,
         const double tolerance = 1e-6,
         const uint _max_iter = 50,
         const miconossmath::normtype normtype = miconossmath::normtype::L2);

  /// method solving a non linear system
  int solve(vecD & lm, matrixD & J, bool & converged, const bool auto_scale_equation = true);

  /// method used to solve a nonlinear system of time depndent equations
  /// allowing to use substepping
  int solveSubstep(vecD & lm,
                   matrixD & J,
                   bool & converged,
                   NLSystemParameters * const sysparams,
                   const std::vector<std::string> & pname,
                   matrixD & Tangent,
                   bool & custom_interruption,
                   double & increment_at_custom_interruption,
                   const uint max_ncut = 1,
                   const bool auto_scale_equation = false,
                   const bool force_substep = false);

protected:
  /// computes the newton update step
  int computeNewtonStep(const vecD & R, const matrixD & J);

  /// custom interrupt while substepping if a given condition is satisfied. Must
  /// override to imlmenet a condition
  virtual int customSubstepInterruption(NLSystemParameters * const sysparams,
                                        bool & custom_interruption_flag);

  /// pointer to the non linear system object
  NLSystem * const _nlsys;
  /// pointer to the non linear system variables
  NLSystemVars * const _sys_vars;
  /// tolerance for convergence
  const double _tolerance;
  /// max number of allowed nonlinear iterations
  const uint _max_iter;
  /// number of equations
  const uint _n_eq;
  /// number of lagrange multipliers
  const uint _n_lm;
  /// total number of unkonws
  const uint _n_total;
  /// the newton increment vector
  vecD _dx;
  /// the selected norm type
  const miconossmath::normtype _normtype;
};
