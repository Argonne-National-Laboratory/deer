#pragma once

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

#include "TNdot.h"
#include "TSdot.h"
#include "aComp.h"
#include "adot.h"
#include "bComp.h"
#include "bdot.h"
#include "equationTests.h"
#include "f_Low.h"
#include "f_aL.h"
#include "f_ab.h"
#include "nlFunBase.h"
#include "prepareVars.h"
#include "qFun.h"
#include "softMax.h"
#include "v1dot.h"
#include "v2dotHigh.h"
#include "v2dotLow.h"
#include "vdot.h"

class GBCavitationNLSystem : public nlFunBase {

public:
  GBCavitationNLSystem(
      const std::string &fname, const Real &D_gb, const Real &h,
      const Real &n_exponent, const Real &E_interface, const Real &E_penalty,
      const Real &G_interface, const Real &eta_sliding,
      const Real &interface_thickness, const Real &beta_exponent,
      const Real &a0, const Real &b0, const Real &b_saturation,
      const Real &sigma_0, const Real &S_thr, const Real &FN,
      const int &vdot_max_type, const int &vdot_type,
      const int &triaxial_vdot_active, const Real &vdot_smooth_factor,
      const bool &cavity_nucleation_on, const bool &cavity_growth_on,
      const bool &triaxial_cavity_growth_on,
      const Real &theta_time_integration);

  // GBCavitationNLSystem(const std::string &f_name, const Real &h, const Real
  // &n, const Real &D);

  DenseVector<Real> computeSystemValue(const io_maps_type &x,
                                       const io_maps_type &params,
                                       const io_maps_type &x_old,
                                       const Real &dt) const;

  GBCavitationNLSystem::io_maps_type
  computeSystemRealValue(const io_maps_type &x_NL_in,
                         const io_maps_type &params, const io_maps_type &x_old,
                         const Real &dt) const;

  DenseMatrix<Real> computeSystemVarJacobian(const io_maps_type &x,
                                             const io_maps_type &params,
                                             const io_maps_type &x_old,
                                             const Real &dt,
                                             const bool &wrt_XNL = true) const;

  std::map<std::string, DenseVector<Real>>
  computeSystemParamGradient(const io_maps_type &x, const io_maps_type &params,
                             const io_maps_type &x_old, const Real &dt) const;

  io_maps_type getRealValueFromNLValue(const io_maps_type &x_NL,
                                       const io_maps_type &x_old) const;

  io_maps_type getNLValueFromRealValue(const io_maps_type &x_Real,
                                       const io_maps_type &x_old) const;

  void returnVolumeRate(const io_maps_type &x, const io_maps_type &params,
                        const io_maps_type &x_old, Real &VL1, Real &VL2,
                        Real &VH1, Real &VH2, Real &Vdot, Real &VLdot,
                        Real &VHdot);
  bool nucleationAboveThreshold(const io_maps_type &params,
                                const io_maps_type &x_old) const {
    return _eqTests.nucleationAboveThreshold(params, x_old);
  };

protected:
  /* INTERFACE NON LINEAR FUNCTIONS */
  equationTests _eqTests;
  f_ab _fab;
  f_aL _faL;
  f_Low _fLow;
  qFun _qLow;
  qFun _qHigh;
  v1dot _v1L;
  v1dot _v1H;
  v2dotLow _v2L;
  v2dotHigh _v2H;
  vdot _vdot;
  TSdot _Ts1dot;
  TSdot _Ts2dot;
  TNdot _TNdot;
  adot _adot;
  bdot _bdot;
  aComp _aComp;
  bComp _bComp;
  prepareVars _prepSNVar;
  const int _vdot_type;
  const int _triaxial_vdot_active;

  void varGradientHelper(const io_maps_type &var_grad,
                         const io_maps_type &myx_dx, DenseMatrix<Real> &J,
                         const unsigned int &idx) const;
};
