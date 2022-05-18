#pragma once
#include "NLParameter.h"
#include "NLPreEquationEvalautionCalc.h"
#include "NLVar.h"
#include "miconosstype.h"

/** Class used to implement a non linear equation i.e. xguess_i -
equation(x1,..xn). Lagrange Multiplier should not be part of the equation as
they will be added automatically. The equation is responsible for computing the
residual and the derivative of the equation w.r.t. variables and parameters.
**/

class Equation
{
public:
  Equation(const unsigned int eq_index,
           NLSystemVars & sysvars,
           NLSystemParameters & sysparams,
           NLPreEquationEvalautionCalc & pre_eval);

  /// the method returning the residual, see getR in RateEquation for an example
  virtual double getR() const = 0;
  /// the method returning the gradient of the equation  w.r.t. variables see
  /// getJRow RateEquation for an example
  virtual vecD getJrow() const = 0;
  /// the method returning the equation index (i.e. the row in the non linear
  /// system)
  uint getIndex() const { return _eq_index; };

  /// method used to check the implemented gradient is correct
  void checkGradient(const double tol = 1e-4, const double eps = 1e-6);

  /// this method must fill the _dequation_dparam vector
  virtual void computeDEquationDP() = 0;

  /// method used to update some constant values that will be used to compute
  /// the equation value
  virtual void updateConstants(){};

  void autoScaleEquation() const;

  ///{@ get and sets methods
  double getDEquationDParam(const uint i) const;
  double getDEquationDParam(const std::string & pname) const;
  void setDEquationDParam(const uint i, const double deqdx);
  void setDEquationDParam(const std::string & pname, const double deqdx);
  ///}@

  /// method to override to set custom scaling rule, default is 1.
  virtual double equationScalingRule() const;

protected:
  /// the equation index in the non linear system
  const uint _eq_index;
  /// the non linear variables
  NLSystemVars & _sys_vars;
  /// the paramters of the system
  NLSystemParameters & _sysparams;
  /// the total number of varaibles (excluding LM)
  const uint _n_vars;
  /// the total number of paramters
  const uint _n_params;
  /// a reference to the method computing shared values
  NLPreEquationEvalautionCalc & _pre_eval;
  /// the gradient of the equation w.r.t. paramters
  vecD _dequation_dparam;
};

/** Class used to implement a rate equation xdot_i =
rate(x1,..xn). The rate and its derivates should be implemented overriding:
computedRate, DComputedRatetDX, and DComputedRatetDP.
The integrated computed values is computed in computedVal which use theta method
to perform time integration. The non linear residual is computed as xguess_i -
computeval. Lagrange Multiplier should not be part of the equation as they will
be added automatically.
**/
class RateEquation : public Equation
{
public:
  RateEquation(const unsigned int eq_index,
               NLSystemVars & sysvars,
               NLSystemParameters & sysparams,
               NLPreEquationEvalautionCalc & pre_eval,
               const double theta = 0);

  /// compute the residual
  double getR() const override final;
  /// compute the dR/dX
  vecD getJrow() const override final;
  /// method to override to update equation constants, in this case we
  /// precalcualte the explicit rate for the theta method
  void updateConstants() override;

protected:
  /// methods to override to implement a rate eqaution
  ///{@
  virtual double computedRate(const bool implicit) const = 0;
  virtual vecD DComputedRatetDx(const bool implicit) const = 0;
  virtual vecD DComputedRatetDP(const bool implicit) const = 0;
  ///@}

private:
  /// computes the derivative of the residaul equation with respect to the
  /// parameters
  void computeDEquationDP() override;
  /// the computed eqaution value after time integration
  double computedVal() const;
  /// the derivative of the computed value
  vecD DComputedValDx() const;
  /// constant used for time integration using the theta method: _theta=0->fully
  /// implicit, _theta=1, fully explicit, _theta=0.5->second order
  const double _theta;
  /// the value of the explicit rate
  double _explicit_rate = 0;
  /// the vector storing dexplicit_rate/dP
  vecD _dexplicit_rate_dp;
};
