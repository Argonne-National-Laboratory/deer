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
        const std::vector<std::string> &value_names, const double n,
        const double h, const double D, const bool use_vl_triax);

  void updateValues(const bool implicit = true) override;
  void updateDerivatives(const bool implicit = true) override;

protected:
  double mFun();
  double alphanFun();
  double betanFun();
  double VLdotFun(const bool implicit);
  vecD dVLdotFundX(const bool implicit);
  double VL1dotFun(const bool implicit);
  vecD dVL1dotFundX(const bool implicit);
  double VL2dotFun(const bool implicit);
  vecD dVL2dotFundX(const bool implicit);
  double VHdotFun(const bool implicit);
  vecD dVHdotFundX(const bool implicit);
  double VH1dotFun(const bool implicit);
  vecD dVH1dotFundX(const bool implicit);
  double VH2dotFun(const bool implicit);
  vecD dVH2dotFundX(const bool implicit);
  double fabFun(const bool implicit);
  vecD dfabFundX(const bool implicit);
  double faLFun(const bool implicit);
  vecD dfaLFundX(const bool implicit);
  double qFun(const bool implicit);
  vecD dqFundX(const bool implicit);
  double qHFun(const bool implicit);
  vecD dqHFundX(const bool implicit);
  double Vdot(const bool implicit);
  vecD dVdotdX(const bool implicit);

  const double _n;
  const double _alpha_n;
  const double _h;
  const double _D;
  const bool _use_vl_triax;
};

class a_res : public RateEquation {
public:
  a_res(const unsigned int eq_index, NLSystemVars &sysvars,
        NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
        const double h, const double a0, const double theta = 0,
        const bool growth_on = true);

  double computedRate(const bool implicit) const override;
  vecD DComputedRatetDx(const bool implicit) const override;
  vecD DComputedRatetDP(const bool implicit) const override;
  // double equationScalingRule() const override;

  const double _h;
  const double _a0;
  const bool _growth_on;
};

class b_res : public RateEquation {
public:
  b_res(const unsigned int eq_index, NLSystemVars &sysvars,
        NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
        const double FN, const double FN_NI, const double S0, const double beta,
        const double b_sat, const double theta = 0,
        const bool nucleation_on = true);

  bool nucleationAboveThreshold(const bool implicit) const;
  bool nucleationIsActive(const bool implicit) const;
  double computedRate(const bool implicit) const override;
  vecD DComputedRatetDx(const bool implicit) const override;
  vecD DComputedRatetDP(const bool implicit) const override;
  // double equationScalingRule() const override;

  const bool _nucleation_on;
  const double _FN;
  const double _FN_NI;
  const double _S0;
  const double _beta;
  const double _b_sat;
};

class TN_res : public RateEquation {
public:
  TN_res(const unsigned int eq_index, NLSystemVars &sysvars,
         NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
         const double thickness, const double E_interface, const double P_mt,
         const double _P_thickness, const double theta = 0);

  double computedRate(const bool implicit) const override;
  vecD DComputedRatetDx(const bool implicit) const override;
  vecD DComputedRatetDP(const bool implicit) const override;
  // double equationScalingRule() const override;

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

  const double _thickness;
  const double _E_interface;
  const double _P_mt;
  const double _P_thickness;
};

class TS_res : public RateEquation {
public:
  TS_res(const uint eq_index, NLSystemVars &sysvars,
         NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
         const uint shear_index, const double thickness,
         const double eta_sliding, const double G_interface,
         const double theta = 0);

public:
  double computedRate(const bool implicit) const override;
  vecD DComputedRatetDx(const bool implicit) const override;
  vecD DComputedRatetDP(const bool implicit) const override;
  // double equationScalingRule() const override;

protected:
  double CS(const bool implicit) const;
  vecD dCSdX(const bool implicit) const;
  double etaFun(const bool implicit) const;
  vecD DetaFunDX(const bool implicit) const;

  const std::string _vname;
  const std::string _udotname;
  const double _thickness;
  const double _eta_sliding;
  const double _G_interface;
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
          const NLSystemParameters &sysparams, const uint n_sys,
          const double a0);

  double gFun() const override;
  vecD dgFun_dx() const override;

  const double _a0;
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
