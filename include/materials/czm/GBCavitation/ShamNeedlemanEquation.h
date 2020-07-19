#pragma once
#include "Equation.h"
#include "InequalityConstraint.h"
#include "NLPreEquationEvalautionCalc.h"
#include "Newton.h"

namespace ShamNeedlemann {

double h_psi(const double psi);

class V_dot : public NLPreEquationEvalautionCalc {

public:
  V_dot(NLSystemVars *const sysvars, const NLSystemParameters *sysparams,
        const std::vector<std::string> &value_names, const bool use_vl_triax);

  void updateValues(const bool implicit = true) override;
  void updateDerivatives(const bool implicit = true) override;

protected:
  double mFun();
  double alphanFun();
  double betanFun();
  double VL2dotFun(const bool implicit);
  vecD dVL2dotFundX(const bool implicit);
  double fabFun(const bool implicit);
  vecD dfabFundX(const bool implicit);
  double faLFun(const bool implicit);
  vecD dfaLFundX(const bool implicit);
  double qFun(const bool implicit);
  vecD dqFundX(const bool implicit);
  double Vdot(const bool implicit);
  vecD dVdotdX(const bool implicit);

  const bool _use_vl_triax;
};

class a_res : public RateEquation {
public:
  a_res(const unsigned int eq_index, NLSystemVars &sysvars,
        NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
        const double theta = 0, const bool growth_on = true);

  double computedRate(const bool implicit) const override;
  vecD DComputedRatetDx(const bool implicit) const override;
  vecD DComputedRatetDP(const bool implicit) const override;
  double equationScalingRule() const override;

  const bool _growth_on;
};

class b_res : public RateEquation {
public:
  b_res(const unsigned int eq_index, NLSystemVars &sysvars,
        NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
        const double theta = 0, const bool nucleation_on = true);

  bool nucleationAboveThreshold(const bool implicit) const;
  bool nucleationIsActive(const bool implicit) const;
  double computedRate(const bool implicit) const override;
  vecD DComputedRatetDx(const bool implicit) const override;
  vecD DComputedRatetDP(const bool implicit) const override;
  double equationScalingRule() const override;

  const bool _nucleation_on;
};

class TN_res : public RateEquation {
public:
  TN_res(const unsigned int eq_index, NLSystemVars &sysvars,
         NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
         const double theta = 0);

  double computedRate(const bool implicit) const override;
  vecD DComputedRatetDx(const bool implicit) const override;
  vecD DComputedRatetDP(const bool implicit) const override;
  double equationScalingRule() const override;

protected:
  double currentJump() const;
  vecD DcurrentJumpDParam() const;
  bool innerPentrationCheck() const;
  double quadraticPenalty() const;
  vecD DQuadraticPenaltyDParam() const;
  double Eeffective() const;
  vecD DEeffectiveDParam() const;
  double CN(const bool implicit) const;
  vecD dCNdX(const bool implicit) const;
  vecD dCNdParam(const bool implicit) const;
};

class TS_res : public RateEquation {
public:
  TS_res(const uint eq_index, NLSystemVars &sysvars,
         NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
         const double theta, const uint shear_index);

public:
  double computedRate(const bool implicit) const override;
  vecD DComputedRatetDx(const bool implicit) const override;
  vecD DComputedRatetDP(const bool implicit) const override;
  double equationScalingRule() const override;

protected:
  double CS(const bool implicit) const;
  vecD dCSdX(const bool implicit) const;
  double etaFun(const bool implicit) const;
  vecD DetaFunDX(const bool implicit) const;

  const std::string _vname;
  const std::string _udotname;
};

class a_lt_b : public InequalityConstraint {
public:
  a_lt_b(const uint lm_index, const NLSystemVars &sys_vars,
         const NLSystemParameters &sysparams, const uint n_sys);

  double gFun() const override;
  vecD dgFun_dx() const override;
};

class a_gt_a0 : public InequalityConstraint {
public:
  a_gt_a0(const uint lm_index, const NLSystemVars &sys_vars,
          const NLSystemParameters &sysparams, const uint n_sys);

  double gFun() const override;
  vecD dgFun_dx() const override;
};

class b_lt_b_old : public InequalityConstraint {
public:
  b_lt_b_old(const uint lm_index, const NLSystemVars &sys_vars,
             const NLSystemParameters &sysparams, const uint n_sys);

  double gFun() const override;
  vecD dgFun_dx() const override;
};

class Solver : public Newton {
public:
  Solver(NLSystem *_nlsys, NLSystemVars *sys_vars, const double tolerance,
         const uint _max_iter, const miconossmath::normtype normtype);

protected:
  int customSubstepInterruption(NLSystemParameters *const sysparams,
                                bool &custom_interruption_flag) override;
};

} // namespace ShamNeedlemann
