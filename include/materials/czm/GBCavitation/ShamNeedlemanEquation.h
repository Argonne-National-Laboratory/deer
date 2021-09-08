//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "Equation.h"
#include "InequalityConstraint.h"
#include "NLPreEquationEvalautionCalc.h"
#include "Newton.h"

namespace ShamNeedlemann {

double h_psi(const double psi);

/// class used to precalcualte the value of Vdot before computing adot and tdot
class V_dot : public NLPreEquationEvalautionCalc {

public:
  V_dot(NLSystemVars *const sysvars, const NLSystemParameters *sysparams,
        const std::vector<std::string> &value_names, const double n,
        const double h, const double D, const bool use_vdot_creep,
        const unsigned int vdot_method, const bool nucleation_on);

  /// update the values and derivatives of Vdot
  void updateValues(const bool implicit = true) override;
  void updateDerivatives(const bool implicit = true) override;

protected:
  /// compute the value of m based on the hydrostatic stress
  double mFun();
  /// compute the alpha as function of the exponent
  double alphanFun();
  /// compute the value of beta
  double betanFun();

  // in what follows X is the vector of nonlinear variables, i.e. a,b,Tn,Ts1,Ts2

  // computes the volume rate as Vdot = max(|VLdot|,|VHdot|)
  double Vdot(const bool implicit);
  // computes dVdot/dX
  vecD dVdotdX(const bool implicit);

  // VOLUME CHANGE FOR LOW TRIAXIALITY REGIME
  /// computes VLdot as VL1dot + VL2dot
  double VLdotFun(const bool implicit);
  /// computes dVLdot dX
  vecD dVLdotFundX(const bool implicit);
  /// computes VL1dot i.e. VL1dot=8 Pi D Tn / q
  double VL1dotFun(const bool implicit);
  /// computes dVl1dot/dX
  vecD dVL1dotFundX(const bool implicit);
  /// computes the volume rate associated to creep
  double VL2dotFun(const bool implicit);
  /// computes dVl1dot/dX
  vecD dVL2dotFundX(const bool implicit);
  /// computes q
  double qFun(const bool implicit);
  /// computes dq/dX
  vecD dqFundX(const bool implicit);

  // VOLUME CHANGE FOR HIGH TRIAXIALITY REGIME
  /// computes VHdot as VH1dot + VH2dot
  double VHdotFun(const bool implicit);
  /// computes dVHdot dX
  vecD dVHdotFundX(const bool implicit);
  /// computes VH1dot i.e. VH1dot=8 Pi D Tn / qH
  double VH1dotFun(const bool implicit);
  /// computes dVH1dot/dX
  vecD dVH1dotFundX(const bool implicit);
  /// computes the volume rate associated to creep
  double VH2dotFun(const bool implicit);
  /// computes dVH2dot/dX
  vecD dVH2dotFundX(const bool implicit);
  /// computes qH
  double qHFun(const bool implicit);
  /// computes dqH/dX
  vecD dqHFundX(const bool implicit);

  /// computes fab=(a/b)^2
  double fabFun(const bool implicit);
  /// computes dfab/dX
  vecD dfabFundX(const bool implicit);
  /// computes faL=(a/(a+L))^2
  double faLFun(const bool implicit);
  /// computes dfaL/dX
  vecD dfaLFundX(const bool implicit);

  /// the value of _n (_n_exponent in the material)
  const double _n;
  ///. the value of alpha
  const double _alpha_n;
  ///. the value of h
  const double _h;
  ///. the diffusivbity value (_D_GB)
  const double _D;
  ///. if true uses the volume increment due to creep (_use_triaxial_growth in
  /// GBCavitation)
  const bool _use_vdot_creep;
  ///. the type of volume increment to use: 0->VLdot, 1->VHdot,
  /// 2->max(|VLdot|,|VHdot|). (_vdot_method in GBCavitaion)
  const unsigned int _vdot_method;
  ///. if true allows nucleation
  const bool _nucleation_on;
};

// in what follows P is the vector of parameters, mostly variables coming from
// the bulk.

/// class implementing the adot equation
class a_res : public RateEquation {
public:
  a_res(const unsigned int eq_index, NLSystemVars &sysvars,
        NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
        const double h, const double a0, const double theta = 0,
        const bool growth_on = true);

  /// computes adot
  double computedRate(const bool implicit) const override;
  /// computes dadot/dX
  vecD DComputedRatetDx(const bool implicit) const override;
  /// computes dadot/dP
  vecD DComputedRatetDP(const bool implicit) const override;

  /// the value of h
  const double _h;
  /// the value of a0
  const double _a0;
  /// if true, cavity growth is active
  const bool _growth_on;
};

/// class implementing the bdot equation
class b_res : public RateEquation {
public:
  b_res(const unsigned int eq_index, NLSystemVars &sysvars,
        NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
        const double FN, const double FN_NI, const double S0, const double beta,
        const double b_sat, const double theta = 0,
        const bool nucleation_on = true);

  /// method checking for the nucleation threshold
  bool nucleationAboveThreshold(const bool implicit) const;
  /// method checking if nucleation is active
  bool nucleationIsActive(const bool implicit) const;
  /// computes bdot
  double computedRate(const bool implicit) const override;
  /// computes dbdot/dX
  vecD DComputedRatetDx(const bool implicit) const override;
  /// computes dbdot/dP
  vecD DComputedRatetDP(const bool implicit) const override;

  /// if true, allows for cavity nucleation
  const bool _nucleation_on;
  /// the value of FN
  const double _FN;
  /// the value of FN_NI
  const double _FN_NI;
  /// the traction normalization value (_sigma_0 in GBCavitation)
  const double _S0;
  /// the cavitation exponent (_beta_exponent in GBCavitation)
  const double _beta;
  /// the saturation value of b
  const double _b_sat;
};

/// class implementing the normal traction equation
class TN_res : public RateEquation {
public:
  TN_res(const unsigned int eq_index, NLSystemVars &sysvars,
         NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
         const double thickness, const double E_interface, const double P_mt,
         const double _P_thickness, const double theta = 0);

  /// computes TN
  double computedRate(const bool implicit) const override;
  /// computes dTn/dX
  vecD DComputedRatetDx(const bool implicit) const override;
  /// computes dTn/dP
  vecD DComputedRatetDP(const bool implicit) const override;

protected:
  /// the value of the normal jump, Un
  double currentJump() const;
  /// dUn/dP
  vecD DcurrentJumpDParam() const;
  /// check for penetration
  bool innerPentrationCheck() const;
  /// computes the quadratic penalty value
  double quadraticPenalty() const;
  /// computes dPenalty/dP
  vecD DQuadraticPenaltyDParam() const;
  /// computes the effective interface modulus, E*Penalty
  double Eeffective() const;
  /// computes dEeffective/dP
  vecD DEeffectiveDParam() const;
  /// computes the stiffness, CN
  double CN(const bool implicit) const;
  /// computes dCN/dX
  vecD dCNdX(const bool implicit) const;
  /// computes dCN/dP
  vecD dCNdParam(const bool implicit) const;

  /// the interface thickness
  const double _thickness;
  /// the interface modulus
  const double _E_interface;
  /// the penalty value when UN=-_P_thickness
  const double _P_mt;
  /// the value of Un at which the penalty is equal to _P_mt
  const double _P_thickness;
};

/// class implementing the shear traction equations
class TS_res : public RateEquation {
public:
  TS_res(const uint eq_index, NLSystemVars &sysvars,
         NLSystemParameters &sysparams, NLPreEquationEvalautionCalc &pre_eval,
         const uint shear_index, const double thickness,
         const double eta_sliding, const double G_interface,
         const double theta = 0);

public:
  /// computes the sliding traction rate dTs
  double computedRate(const bool implicit) const override;
  /// computes dTs/dX
  vecD DComputedRatetDx(const bool implicit) const override;
  /// computes dTs/dP
  vecD DComputedRatetDP(const bool implicit) const override;

protected:
  /// computes the shear stiffness CS
  double CS(const bool implicit) const;
  /// computes dCS/dP
  vecD dCSdX(const bool implicit) const;
  /// computes the sliding viscosity modifier eta
  double etaFun(const bool implicit) const;
  /// computesthe deta/dX
  vecD DetaFunDX(const bool implicit) const;

  /// the name of the shear traction: Ts1, or Ts2
  const std::string _vname;
  /// the name of the dispalcement jump rate, udot_S1, or udot_2
  const std::string _udotname;
  /// the interface thickness
  const double _thickness;
  /// the initial sliding viscosity
  const double _eta_sliding;
  /// the interface shear modulus
  const double _G_interface;
};

// class implementing the inequality constraint a<b
class a_lt_b : public InequalityConstraint {
public:
  a_lt_b(const uint lm_index, const NLSystemVars &sys_vars,
         const NLSystemParameters &sysparams, const uint n_sys);

  /// computes g=a/b-ab_max
  double gFun() const override;
  /// computes dg/dX
  vecD dgFun_dx() const override;
};

// class implementing the inequality constraint a>=a0
class a_gt_a0 : public InequalityConstraint {
public:
  a_gt_a0(const uint lm_index, const NLSystemVars &sys_vars,
          const NLSystemParameters &sysparams, const uint n_sys,
          const double a0);

  /// computes a0/a-1
  double gFun() const override;
  /// computes dg/dX
  vecD dgFun_dx() const override;
  /// the value of a0
  const double _a0;
};

// class implementing the inequality constraint b>=b_old
class b_lt_b_old : public InequalityConstraint {
public:
  b_lt_b_old(const uint lm_index, const NLSystemVars &sys_vars,
             const NLSystemParameters &sysparams, const uint n_sys);

  /// computes g=b-bold
  double gFun() const override;
  /// computes dg/dx
  vecD dgFun_dx() const override;
};

/// custom newton solver allowing for custom substep interruption
class Solver : public Newton {
public:
  Solver(NLSystem *_nlsys, NLSystemVars *sys_vars, const double tolerance,
         const uint _max_iter, const miconossmath::normtype normtype);

protected:
  /// method checking if we need to interrupt while substepping.
  int customSubstepInterruption(NLSystemParameters *const sysparams,
                                bool &custom_interruption_flag) override;
};

} // namespace ShamNeedlemann
