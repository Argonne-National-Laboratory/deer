#pragma once
#include "Equation.h"
#include "InequalityConstraint.h"
#include "NLPreEquationEvalautionCalc.h"
#include "Newton.h"

namespace ShamNeedlemann {

double h_psi(const double psi) {
  const double c = std::cos(psi);
  const double s = std::sin(psi);
  return (1. / (1 - c) - c / 2.) * 1 / s;
}
} // namespace ShamNeedlemann

class V_dot : public NLPreEquationEvalautionCalc {
public:
  using NLPreEquationEvalautionCalc::NLPreEquationEvalautionCalc;

  double mFun() { return _sysparams->getValue("Sh") >= 0 ? 1 : -1; }

  double alphanFun() { return 3. / (2. * _sysparams->getValue("n")); }

  double betanFun() {
    const double g1 = std::log(3) - 2. / 3.;
    const double gm1 = 2. * M_PI / (9. * std::sqrt(3));
    const double g = _sysparams->getValue("Sh") >= 0 ? g1 : gm1;
    const double n = _sysparams->getValue("n");
    return (n - 1) * (n + g) / (n * n);
  }

  double VL2dotFun(const bool implicit) {
    const double svm = _sysparams->getValue("Svm");
    const double sh = _sysparams->getValue("Sh");
    double vL2dot = 0;
    if (svm != 0 && sh != 0) {
      const double m = mFun();
      const double triax = sh / svm;
      const double a = _sys_vars->getValueImplicit("a", implicit);
      const double alpha_n = alphanFun();
      const double betan = betanFun();
      const double n = _sysparams->getValue("n");

      double vL2dot = 2 * _sysparams->getValue("edot") * (a * a * a) * M_PI *
                      _sysparams->getValue("h");

      if (std::abs(triax) >= 1)
        vL2dot *= m * std::pow(alpha_n * std::abs(triax) + betan, n);
      else
        vL2dot *= std::pow(alpha_n + betan, n) * triax;
    }
    return vL2dot;
  }

  vecD dVL2dotFundX(const bool implicit) {
    const double svm = _sysparams->getValue("Svm");
    const double sh = _sysparams->getValue("Sh");
    vecD dvL2dotdx(_n_vars);
    if (svm != 0 && sh != 0) {
      const double m = mFun();
      const double triax = sh / svm;
      const double a = _sys_vars->getValueImplicit("a", implicit);
      const double alpha_n = alphanFun();
      const double betan = betanFun();
      const double n = _sysparams->getValue("n");

      dvL2dotdx[0] = 6 * _sysparams->getValue("edot") * (a * a) * M_PI *
                     _sysparams->getValue("h");

      if (std::abs(triax) >= 1)
        dvL2dotdx[0] *= m * std::pow(alpha_n * std::abs(triax) + betan, n);
      else
        dvL2dotdx[0] *= std::pow(alpha_n + betan, n) * triax;
    }
    return dvL2dotdx;
  }

  double fabFun(const bool implicit) {
    const double a = _sys_vars->getValueImplicit("a", implicit);
    const double b = _sys_vars->getValueImplicit("b", implicit);
    return a * a / (b * b);
  }

  vecD dfabFundX(const bool implicit) {
    vecD dfabdx(_n_vars);
    const double a = _sys_vars->getValueImplicit("a", implicit);
    const double b = _sys_vars->getValueImplicit("b", implicit);

    dfabdx[0] = 2 * a / (b * b);
    dfabdx[1] = -2 * a * a / (b * b * b);

    return dfabdx;
  }

  double faLFun(const bool implicit) {

    const double edot = _sysparams->getValue("edot");
    const double svm = _sysparams->getValue("Svm");
    double L = 0;
    double FL = 0;

    if (edot != 0 && svm != 0) {
      const double a = _sys_vars->getValueImplicit("a", implicit);
      L = std::pow(_sysparams->getValue("D") * svm / edot, 1. / 3.);
      FL = a * a / ((a + 1.5 * L) * (a + 1.5 * L));
    }
    return FL;
  }

  vecD dfaLFundX(const bool implicit) {
    vecD dfaldx(_n_vars);
    const double edot = _sysparams->getValue("edot");
    const double svm = _sysparams->getValue("Svm");
    double L = 0;
    if (edot != 0 && svm != 0) {

      const double a = _sys_vars->getValueImplicit("a", implicit);
      L = std::pow(_sysparams->getValue("D") * svm / edot, 1. / 3.);
      dfaldx[0] = 24. * L * a / std::pow(2 * a + 3 * L, 3);
    }
    return dfaldx;
  }

  double qFun(const bool implicit) {
    // const double p = 30;
    // const double f = std::pow(
    //     std::pow(fabFun(implicit), p) + std::pow(faLFun(implicit), p), 1. /
    //     p);
    double f = 0;
    const double fab = fabFun(implicit);
    const double faL = fabFun(implicit);
    if (fab >= faL)
      f = fab;
    else
      f = faL;

    return -2. * std::log(f) - (1. - f) * (3. - f);
  }

  vecD dqFundX(const bool implicit) {

    const double fab = fabFun(implicit);
    const double faL = faLFun(implicit);
    vecD dqdx(_n_vars);

    // const double p = 30;
    // const double fabP = std::pow(fab, p);
    // const double faLP = std::pow(faL, p);
    // const double f = std::pow(fabP + faLP, 1. / p);
    //
    // const double prefactor = std::pow(fabP + faLP, (1. - p) / p);
    //
    // vecD dfabdX = dfabFundX(implicit);
    // vecD dfaLdX = dfaLFundX(implicit);
    //
    //
    // for (uint i = 0; i < _n_vars; i++)
    //   dqdx[i] = prefactor *
    //             (fabP / fab * dfabdX[i] + std::pow(faL, p - 1.) * dfaLdX[i]);

    double f;
    if (fab >= faL) {
      f = fab;
      dqdx = dfabFundX(implicit);
    } else {
      f = faL;
      dqdx = dfaLFundX(implicit);
    }
    const double dqdf = -2. * f - 2. / f + 4.;
    for (uint i = 0; i < _n_vars; i++)
      dqdx[i] *= dqdf;

    return dqdx;
  }

  double Vdot(const bool implicit) {
    double vdot = 8 * M_PI * _sysparams->getValue("D") *
                  _sys_vars->getValueImplicit("Tn", implicit) / qFun(implicit);

    return vdot + VL2dotFun(implicit);
  }

  vecD dVdotdX(const bool implicit) {
    vecD dvdot_dx(_n_vars);
    const double num = _sys_vars->getValueImplicit("Tn", implicit);
    const double prefactor = 8 * M_PI * _sysparams->getValue("D");

    vecD dnum_dx(_n_vars);
    dnum_dx[2] = 1;

    const double q = qFun(implicit);
    const vecD dqdx = dqFundX(implicit);
    const vecD dvL2dx = dVL2dotFundX(implicit);

    for (uint i = 0; i < _n_vars; i++)
      dvdot_dx[i] =
          prefactor * (dnum_dx[i] / q - num * dqdx[i] / (q * q)) + dvL2dx[i];

    return dvdot_dx;
  }

  void updateValues(const bool implicit = true) override {
    setValue("vdot", Vdot(implicit), implicit);
  }

  void updateDerivatives(const bool implicit = true) override {
    setDValueDX("vdot", dVdotdX(implicit), implicit);
  }
};

class a_res : public RateEquation {
public:
  a_res(const unsigned int eq_index, NLSystemVars &sysvars,
        NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
        const double theta = 0, const bool growth_on = true)
      : RateEquation(eq_index, sysvars, sysparams, pre_eval, theta),
        _growth_on(growth_on) {}

  const bool _growth_on;

  double computedRate(const bool implicit) const override {

    // a_dot = V_dot/(4*pi*a^2*h)
    double a_dot = 0;

    if (_growth_on) {
      const double udot = _sysparams.getValue("udot_N");
      if (udot > 0 || (udot < 0 && _sys_vars.getValueOld("a") >
                                       _sysparams.getValue("a0"))) {
        const double a = _sys_vars.getValueImplicit("a", implicit);
        a_dot = _pre_eval.getValue("vdot", implicit) /
                (4. * M_PI * _sysparams.getValue("h") * a * a);
      }
    }
    return a_dot;
  }

  vecD DComputedRatetDx(const bool implicit) const override {
    vecD dadot_dx(_n_vars);
    if (_growth_on) {
      const double udot = _sysparams.getValue("udot_N");
      if (udot > 0 || (udot < 0 && _sys_vars.getValueOld("a") >
                                       _sysparams.getValue("a0"))) {
        const double a = _sys_vars.getValueImplicit("a", implicit);
        const double vdot = _pre_eval.getValue("vdot", implicit);
        const vecD dvdot_dx = _pre_eval.getDValueDX("vdot", implicit);
        const double h = _sysparams.getValue("h");

        double g = (4. * M_PI * h * a * a);
        vecD dg_dx(_n_vars);
        dg_dx[0] = (8. * M_PI * h * a);

        for (uint i = 0; i < _n_vars; i++)
          dadot_dx[i] = dvdot_dx[i] / g - vdot * dg_dx[i] / (g * g);
      }
    }
    return dadot_dx;
  }

  vecD DComputedRatetDP(const bool implicit) const override {
    return vecD(_n_params);
  }
};

class b_res : public RateEquation {
public:
  b_res(const unsigned int eq_index, NLSystemVars &sysvars,
        NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
        const double theta = 0, const bool nucleation_on = true)
      : RateEquation(eq_index, sysvars, sysparams, pre_eval, theta),
        _nucleation_on(nucleation_on) {}

  const bool _nucleation_on;

  bool nucleationAboveThreshold(const bool implicit) const {
    bool active = false;
    const double tn = _sys_vars.getValueImplicit("Tn", implicit);
    if (tn > 0)
      active = (std::pow(tn / _sysparams.getValue("S0"),
                         _sysparams.getValue("beta")) *
                _sysparams.getValue("e")) > (1. / _sysparams.getValue("FN_NI"));
    return active;
  }

  bool nucleationIsActive(const bool implicit) const {
    bool active = _sysparams.getValue("nucleation_is_active");
    if (active) {
      active =
          active && (_sys_vars.getValueOld("b") > _sysparams.getValue("b_sat"));
      active = active && (_sys_vars.getValueImplicit("Tn", implicit) > 0);
    } else
      active = nucleationAboveThreshold(implicit);
    return active;
  }

  double computedRate(const bool implicit) const override {

    double bdot = 0;
    if (_nucleation_on && nucleationIsActive(implicit)) {
      const double b = _sys_vars.getValueImplicit("b", implicit);
      bdot = -M_PI * (b * b * b) * _sysparams.getValue("FN") *
             std::pow(_sys_vars.getValueImplicit("Tn", implicit) /
                          _sysparams.getValue("S0"),
                      _sysparams.getValue("beta")) *
             _sysparams.getValue("edot");
    }
    return bdot;
  }

  vecD DComputedRatetDx(const bool implicit) const override {
    vecD dbdot_dx(_n_vars);
    if (_nucleation_on && nucleationIsActive(implicit)) {
      const double b = _sys_vars.getValueImplicit("b", implicit);
      const double T = _sys_vars.getValueImplicit("Tn", implicit);
      const double S0 = _sysparams.getValue("S0");
      const double beta = _sysparams.getValue("beta");
      const double edot = _sysparams.getValue("edot");
      const double FN = _sysparams.getValue("FN");
      dbdot_dx[1] = -3. * M_PI * (b * b) * FN * std::pow(T / S0, beta) * edot;
      dbdot_dx[2] = -beta * M_PI * (b * b * b) * FN *
                    std::pow(T / S0, beta - 1.) / S0 * edot;
    }
    return dbdot_dx;
  }

  vecD DComputedRatetDP(const bool implicit) const override {
    return vecD(_n_params);
  }
};

class TN_res : public RateEquation {
public:
  using RateEquation::RateEquation;

  double CN(const bool implicit) const {

    return _sysparams.getValue("thickness") /
           (_sysparams.getValue("E") *
            (1. - _sys_vars.getValueImplicit("a", implicit) /
                      _sys_vars.getValueImplicit("b", implicit)));
  }

  vecD dCNdX(const bool implicit) const {
    vecD dCN_dx(_n_vars);
    const double a = _sys_vars.getValueImplicit("a", implicit);
    const double b = _sys_vars.getValueImplicit("b", implicit);
    const double thickness = _sysparams.getValue("thickness");
    const double E = _sysparams.getValue("E");

    const double prefactor = thickness / (E * (b - a) * (b - a));

    dCN_dx[0] = b * prefactor;
    dCN_dx[1] = -a * prefactor;

    return dCN_dx;
  }

  double computedRate(const bool implicit) const override {
    const double b = _sys_vars.getValueImplicit("b", implicit);
    double Tdot = (_sysparams.getValue("udot_N") -
                   _pre_eval.getValue("vdot", implicit) / (M_PI * b * b)) /
                  CN(implicit);

    return Tdot;
  }

  vecD DComputedRatetDx(const bool implicit) const override {
    vecD dTdot_dx(_n_vars);
    const double b = _sys_vars.getValueImplicit("b", implicit);
    const double u_dot = _sysparams.getValue("udot_N");
    const double vdot = _pre_eval.getValue("vdot", implicit);
    const vecD dvdot_dx = _pre_eval.getDValueDX("vdot", implicit);
    const double cn = CN(implicit);
    const vecD dcn_dx = dCNdX(implicit);

    double g = M_PI * b * b;
    vecD dg_dx(_n_vars);
    dg_dx[1] = M_PI * 2 * b;

    const double term1 = (u_dot - vdot / g);
    vecD dterm1_dx(_n_vars);
    for (uint i = 0; i < _n_vars; i++)
      dterm1_dx[i] = -dvdot_dx[i] / g + vdot * dg_dx[i] / (g * g);

    for (uint i = 0; i < _n_vars; i++)
      dTdot_dx[i] = dterm1_dx[i] / cn - term1 * dcn_dx[i] / (cn * cn);

    return dTdot_dx;
  }

  vecD DComputedRatetDP(const bool implicit) const override {
    vecD deq_dparam(_n_params);
    deq_dparam[_sysparams.getParamIndex("udot_N")] = 1. / CN(implicit);
    return deq_dparam;
  }

  double equationScalingRule() const override {
    return std::max(std::abs(_sys_vars.getValueOld(_eq_index)), 10.);
  }
};

class TS_res : public RateEquation {
public:
  TS_res(const uint eq_index, NLSystemVars &sysvars,
         NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
         const double theta, const uint shear_index)
      : RateEquation(eq_index, sysvars, sysparams, pre_eval, theta),
        _vname("Ts" + std::to_string(shear_index)),
        _udotname("udot_S" + std::to_string(shear_index)) {}

  const std::string _vname;
  const std::string _udotname;

  double CS(const bool implicit) const {

    return _sysparams.getValue("thickness") /
           (_sysparams.getValue("G") *
            (1. - _sys_vars.getValueImplicit("a", implicit) /
                      _sys_vars.getValueImplicit("b", implicit)));
  }

  vecD dCSdX(const bool implicit) const {
    vecD dCS_dx(_n_vars);
    const double a = _sys_vars.getValueImplicit("a", implicit);
    const double b = _sys_vars.getValueImplicit("b", implicit);
    const double thickness = _sysparams.getValue("thickness");
    const double G = _sysparams.getValue("G");

    const double prefactor = thickness / (G * (b - a) * (b - a));

    dCS_dx[0] = b * prefactor;
    dCS_dx[1] = -a * prefactor;

    return dCS_dx;
  }

  double etaFun(const bool implicit) const {
    double eta = _sysparams.getValue("eta_sliding");
    const double a_b = _sys_vars.getValueImplicit("a", implicit) /
                       _sys_vars.getValueImplicit("b", implicit);

    if (a_b > 0.5)
      eta *= 2. * (-a_b + 1.);
    return eta;
  }

  vecD DetaFunDX(const bool implicit) const {
    vecD deta_dx(_n_vars);

    const double a_b = _sys_vars.getValueImplicit("a", implicit) /
                       _sys_vars.getValueImplicit("b", implicit);

    if (a_b > 0.5) {
      const double b = _sys_vars.getValueImplicit("b", implicit);
      const double temp = _sysparams.getValue("eta_sliding") * 2.;
      deta_dx[0] = -temp / b;
      deta_dx[1] = temp * a_b / b;
    }
    return deta_dx;
  }

  double computedRate(const bool implicit) const override {
    double Tdot =
        (_sysparams.getValue(_udotname) -
         _sys_vars.getValueImplicit(_vname, implicit) / etaFun(implicit)) /
        CS(implicit);

    return Tdot;
  }

  vecD DComputedRatetDx(const bool implicit) const override {
    vecD dTdot_dx(_n_vars);

    const double u_dot = _sysparams.getValue(_udotname);
    const double ts = _sys_vars.getValueImplicit(_vname, implicit);
    vecD dts_dx(_n_vars);
    dts_dx[_eq_index] = 1;

    const double cs = CS(implicit);
    const vecD dcs_dx = dCSdX(implicit);

    const double g = etaFun(implicit);
    const vecD dg_dx = DetaFunDX(implicit);

    const double term1 = (u_dot - ts / g);
    vecD term1_dx(_n_vars);
    for (uint i = 0; i < _n_vars; i++)
      term1_dx[i] = -dts_dx[i] / g + ts * dg_dx[i] / (g * g);

    for (uint i = 0; i < _n_vars; i++)
      dTdot_dx[i] = -term1 / (cs * cs) * dcs_dx[i] + term1_dx[i] / cs;

    return dTdot_dx;
  }

  vecD DComputedRatetDP(const bool implicit) const override {
    vecD deq_dparam(_n_params);
    deq_dparam[_sysparams.getParamIndex(_udotname)] = 1. / CS(implicit);
    return deq_dparam;
  }

  double equationScalingRule() const override {
    return std::max(std::abs(_sys_vars.getValueOld(_vname)), 1.);
  }
};

class a_lt_b : public InequalityConstraint {
public:
  using InequalityConstraint::InequalityConstraint;
  double gFun() const override {
    return _sys_vars.getValue("a") / _sys_vars.getValue("b") -
           _sysparams.getValue("max_ab");
  }

  vecD dgFun_dx() const override {
    vecD Rgrad(_n_vars, 0);
    double a = _sys_vars.getValue("a");
    double da_dascaled = _sys_vars.getDVarDVarScaled("a");
    double b = _sys_vars.getValue("b");
    double db_dbscaled = _sys_vars.getDVarDVarScaled("b");
    Rgrad[0] = da_dascaled / b;
    Rgrad[1] = -a / (b * b) * db_dbscaled;
    return Rgrad;
  }
};

class a_gt_a0 : public InequalityConstraint {
public:
  using InequalityConstraint::InequalityConstraint;

  double gFun() const override {
    return _sysparams.getValue("a0") / _sys_vars.getScalingFactor("a") -
           _sys_vars.getValueScaled("a");
  }

  vecD dgFun_dx() const override {
    vecD Rgrad(_n_vars, 0);
    Rgrad[0] = -1;
    return Rgrad;
  }
};

class b_lt_b_old : public InequalityConstraint {
public:
  using InequalityConstraint::InequalityConstraint;

  double gFun() const override {
    return _sys_vars.getValueScaled("b") - _sys_vars.getValueOldScaled("b");
  }

  vecD dgFun_dx() const override {
    vecD Rgrad(_n_vars, 0);
    Rgrad[1] = 1;
    return Rgrad;
  }
};

class ShamNeedlemannNewtonSolver : public Newton {
public:
  using Newton::Newton;

protected:
  bool customSubstepInterruption(NLSystemParameters *const sysparams) override {
    const double D = _sys_vars->getValue("a") / _sys_vars->getValue("b");
    const double D_old =
        _sys_vars->getValueOld("a") / _sys_vars->getValueOld("b");
    const double D_rate = (D - D_old) / sysparams->getValue("dt");
    const double residual_life = D_rate > 0 ? (1. - D) / D_rate : 1e6;

    sysparams->setValue("residual_life", residual_life);

    if (D > sysparams->getValue("max_damage")) {
      sysparams->setValue("element_failed", 1.);
      return true;
    }

    if (D_rate > 0)
      if (residual_life <
          sysparams->getValue("minimum_allowed_residual_life")) {
        sysparams->setValue("element_failed", 1.);
        return true;
      }

    if (_sys_vars->getValue("Tn") >
        sysparams->getValue("maximum_allowed_opening_traction")) {
      sysparams->setValue("element_failed", 1.);
      return true;
    }

    sysparams->setValue("element_failed", 0.);
    return false;
  }
};
