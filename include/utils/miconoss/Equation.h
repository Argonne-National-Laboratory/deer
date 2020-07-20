#pragma once
#include "NLParameter.h"
#include "NLPreEquationEvalautionCalc.h"
#include "NLVar.h"
#include "miconosstype.h"

class Equation {
public:
  Equation(const unsigned int eq_index, NLSystemVars &sysvars,
           NLSystemParameters &sysparams,
           NLPreEquationEvalautionCalc &pre_eval);
  virtual double getR() const = 0;
  virtual vecD getJrow() const = 0;
  uint getIndex() const { return _eq_index; };

  void checkGradinet(double tol = 1e-4, double eps = 1e-6);

  // this function must fill the _dequation_dparam vector
  virtual void computeDEquationDP() = 0;
  virtual void updateConstants(){};
  void autoScaleEquation() const;

  ///{@ get and sets methods
  double getDEquationDParam(const uint i) const;
  double getDEquationDParam(const std::string &pname) const;
  void setDEquationDParam(const uint i, const double deqdx);
  void setDEquationDParam(const std::string &pname, const double deqdx);
  ///}@

  /// method to override to set custom scaling rule, default is 1.
  virtual double equationScalingRule() const;

protected:
  const uint _eq_index;
  NLSystemVars &_sys_vars;
  NLSystemParameters &_sysparams;
  const uint _n_vars;
  const uint _n_params;
  NLPreEquationEvalautionCalc &_pre_eval;
  vecD _dequation_dparam;
};

class RateEquation : public Equation {
public:
  RateEquation(const unsigned int eq_index, NLSystemVars &sysvars,
               NLSystemParameters &sysparams,
               NLPreEquationEvalautionCalc &pre_eval, const double theta = 0);

  double getR() const override final;
  vecD getJrow() const override final;

  void updateConstants() override;

protected:
  virtual double computedRate(const bool implicit) const = 0;
  virtual vecD DComputedRatetDx(const bool implicit) const = 0;
  virtual vecD DComputedRatetDP(const bool implicit) const = 0;

private:
  void computeDEquationDP() override;
  double computedVal() const;
  vecD DComputedValDx() const;
  const double _theta;
  double _explicit_rate = 0;
  vecD _dexplicit_rate_dp;
};
