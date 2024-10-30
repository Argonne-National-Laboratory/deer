#pragma once

#include "CZMComputeLocalTractionTotalBase.h"
#include "RankTwoTensor.h"

namespace MathUtils
{
template <typename T>
T
heaviside(T x)
{
  return x > T(0) ? T(1) : T(0);
}

template <typename T>
T
macaulay(T x)
{
  return heaviside(x) * x;
}
} // namespace MathUtils

/**
 * Implementation of the grain boundary cavitation model
 * Nassif, Omar, et al. "Combined crystal plasticity and grain boundary modeling of creep in
 * ferritic-martensitic steels: I. Theory and implementation." Modelling and Simulation in Materials
 * Science and Engineering 27.7 (2019): 075009.
 **/
class GrainBoundaryCavitation : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  GrainBoundaryCavitation(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();

  virtual void computeQpProperties();

  virtual void computeDamage();

  virtual void computeInterfaceTractionAndDerivatives();

  virtual void computeDamageDrivingForces();

  virtual void computeTimeStepLimit();

  /// @{ Bulk properties
  const MaterialProperty<RankTwoTensor> & _stress_primary;
  const MaterialProperty<RankTwoTensor> & _stress_secondary;
  const MaterialProperty<RankTwoTensor> & _creep_strain_primary;
  const MaterialProperty<RankTwoTensor> & _creep_strain_secondary;
  MaterialProperty<Real> & _ec;
  const MaterialProperty<Real> & _ec_old;
  MaterialProperty<Real> & _delta_ec;
  MaterialProperty<Real> & _sigma_vm;
  MaterialProperty<Real> & _sigma_h;
  MaterialProperty<Real> & _Tn;
  MaterialProperty<Real> & _eta;
  /// @}

  /// @{ Damage driving forces
  MaterialProperty<Real> & _delta_N;
  const MaterialProperty<Real> & _delta_N_old;
  MaterialProperty<Real> & _delta_V;
  const MaterialProperty<Real> & _delta_V_old;
  ///@}

  /// @{ Internal state variables
  MaterialProperty<Real> & _a;
  const MaterialProperty<Real> & _a_old;
  MaterialProperty<Real> & _b;
  const MaterialProperty<Real> & _b_old;
  MaterialProperty<Real> & _D;
  const MaterialProperty<Real> & _D_old;
  ///@}

  /// @{ Initial conditions
  const MaterialProperty<Real> & _a0;
  const MaterialProperty<Real> & _a0_old;
  const MaterialProperty<Real> & _b0;
  const MaterialProperty<Real> & _b0_old;
  /// @}

  /// @{ Model parameters. See MooseDocs for their descriptions.
  const Real _psi;
  const MaterialProperty<Real> & _D_GB;
  const Real _n;
  const MaterialProperty<Real> & _E;
  const MaterialProperty<Real> & _G;
  const MaterialProperty<Real> & _w;
  const MaterialProperty<Real> & _eta_s;
  const Real _p;
  const Real _eps;
  const Real _P;
  const Real _gamma;
  const MaterialProperty<Real> & _T0;
  const MaterialProperty<Real> & _FN;
  const MaterialProperty<Real> & _Nc;
  MaterialProperty<Real> & _nucleating;
  const MaterialProperty<Real> & _nucleating_old;
  /// @}

  MaterialProperty<Real> & _dt_max;

  const Real _delta_D_max;
  const Real _dt_cutback_factor;

  const bool _diffusion_growth;
  const bool _creep_growth;

  const bool _fixed_triaxiality_state;
  const MooseEnum _triaxiality_state;

private:
  Real cavityVolume() const;
  Real cavityShapeFactor() const;
  Real degradation() const;
  Real effectiveDamage() const;
  Real normalSeparation() const;
  std::pair<Real, Real> penetrationPenalty(Real ju_n) const;
  Real normalStiffness() const;
  Real tangentialStiffness() const;
  Real viscousStiffness() const;
};
