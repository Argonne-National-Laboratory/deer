#include "GrainBoundaryCavitation.h"
#include "MathUtils.h"

registerMooseObject("DeerApp", GrainBoundaryCavitation);

InputParameters
GrainBoundaryCavitation::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();

  params.addClassDescription("Grain boundary cavitation model");

  params.addRequiredParam<MaterialPropertyName>("a0", "The initial average cavity radius");
  params.addRequiredParam<MaterialPropertyName>("b0", "The initial average cavity half spacing");

  params.addRequiredParam<Real>("psi",
                                "The angle (in degrees) between the tangent of the tip of the "
                                "spherical-cap shaped void and the (local) horizontal plane");
  params.addParam<MaterialPropertyName>("D_GB", "D_GB", "Grain boundary diffusivity");
  params.addRequiredParam<Real>("n", "Bulk power law creep exponent");
  params.addParam<MaterialPropertyName>("E", "E", "Grain boundary Young's modulus");
  params.addParam<MaterialPropertyName>("G", "G", "Grain boundary shear modulus");
  params.addParam<MaterialPropertyName>("w", "w", "Grain boundary thickness");
  params.addRequiredParam<Real>("eps", "Residual stiffness");
  params.addRequiredParam<Real>("P", "Penalty to prevent interpenatration");
  params.addRequiredParam<Real>("gamma", "Cavity nucleation exponent");
  params.addParam<MaterialPropertyName>(
      "T0", "T0", "Reference traction (used to normalize the interface traction)");
  params.addParam<MaterialPropertyName>(
      "FN", "FN", "Grain boundary reference cavity nucleation rate");
  params.addParam<MaterialPropertyName>(
      "Nc", "Nc", "Threshold number density for cavity nucleation");
  params.addParam<MaterialPropertyName>("eta_s", "eta_s", "GB sliding viscosity");
  params.addParam<Real>("p", 2, "The exponent for degradation as a function of damage");
  params.addParam<Real>("delta_D_max",
                        0.05,
                        "Maximum damage increment per step. If the damage increment exceeds this "
                        "threshold, the material timestep limit will be reduced accordingly.");
  params.addParam<Real>("timestep_cutback_factor",
                        0.1,
                        "Factor used to reduce the timestep size when damage increases too fast.");
  params.addParam<bool>(
      "growth_due_to_diffusion", true, "Whether to account for void growth due to diffusion");
  params.addParam<bool>(
      "growth_due_to_creep", true, "Whether to account for void growth due to creep");

  MooseEnum triaxiality_state("LOW MEDIUM HIGH", "LOW");
  params.addParam<MooseEnum>("fixed_triaxiality",
                             triaxiality_state,
                             "Fix the stress triaxiality eta to improve convergence. LOW for eta < "
                             "0, MEDIUM for 0 <= eta < 1, HIGH for eta > 1");

  return params;
}

GrainBoundaryCavitation::GrainBoundaryCavitation(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    // dependent bulk properties
    _stress_primary(getMaterialProperty<RankTwoTensor>(_base_name + "cauchy_stress")),
    _stress_secondary(getNeighborMaterialProperty<RankTwoTensor>(_base_name + "cauchy_stress")),
    _creep_strain_primary(getMaterialProperty<RankTwoTensor>(_base_name + "inelastic_strain")),
    _creep_strain_secondary(
        getNeighborMaterialProperty<RankTwoTensor>(_base_name + "inelastic_strain")),
    _ec(declareProperty<Real>(_base_name + "equivalent_creep_strain")),
    _ec_old(getMaterialPropertyOldByName<Real>(_base_name + "equivalent_creep_strain")),
    _delta_ec(declareProperty<Real>(_base_name + "equivalent_creep_strain_increment")),
    _sigma_vm(declareProperty<Real>(_base_name + "vonmises_stress")),
    _sigma_h(declareProperty<Real>(_base_name + "hydrostatic_stress")),
    _Tn(declareProperty<Real>(_base_name + "normal_traction")),
    _eta(declareProperty<Real>(_base_name + "stress_triaxiality")),
    // damage driving forces
    _delta_N(declareProperty<Real>(_base_name + "cavity_number_density_inc")),
    _delta_N_old(getMaterialPropertyOldByName<Real>(_base_name + "cavity_number_density_inc")),
    _delta_V(declareProperty<Real>(_base_name + "cavity_volume_inc")),
    _delta_V_old(getMaterialPropertyOldByName<Real>(_base_name + "cavity_volume_inc")),
    // internal state variables
    _a(declareProperty<Real>(_base_name + "average_cavity_radius")),
    _a_old(getMaterialPropertyOldByName<Real>(_base_name + "average_cavity_radius")),
    _b(declareProperty<Real>(_base_name + "average_cavity_half_spacing")),
    _b_old(getMaterialPropertyOldByName<Real>(_base_name + "average_cavity_half_spacing")),
    _D(declareProperty<Real>(_base_name + "damage")),
    _D_old(getMaterialPropertyOldByName<Real>(_base_name + "damage")),
    // initial conditions
    _a0(getMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("a0"))),
    _a0_old(getMaterialPropertyOldByName<Real>(_base_name + getParam<MaterialPropertyName>("a0"))),
    _b0(getMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("b0"))),
    _b0_old(getMaterialPropertyOldByName<Real>(_base_name + getParam<MaterialPropertyName>("b0"))),
    // model parameters
    _psi(getParam<Real>("psi") / 180 * M_PI),
    _D_GB(getMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("D_GB"))),
    _n(getParam<Real>("n")),
    _E(getMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("E"))),
    _G(getMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("G"))),
    _w(getMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("w"))),
    _eta_s(getMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("eta_s"))),
    _p(getParam<Real>("p")),
    _eps(getParam<Real>("eps")),
    _P(getParam<Real>("P")),
    _gamma(getParam<Real>("gamma")),
    _T0(getMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("T0"))),
    _FN(getMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("FN"))),
    _Nc(getMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("Nc"))),
    _nucleating(declareProperty<Real>(_base_name + "nucleating")),
    _nucleating_old(getMaterialPropertyOldByName<Real>(_base_name + "nucleating")),
    _dt_max(declareProperty<Real>(_base_name + "material_timestep_limit")),
    _delta_D_max(getParam<Real>("delta_D_max")),
    _dt_cutback_factor(getParam<Real>("timestep_cutback_factor")),
    _diffusion_growth(getParam<bool>("growth_due_to_diffusion")),
    _creep_growth(getParam<bool>("growth_due_to_creep")),
    // should we fix stress triaxiality state?
    _fixed_triaxiality_state(isParamSetByUser("fixed_triaxiality")),
    _triaxiality_state(getParam<MooseEnum>("fixed_triaxiality"))
{
}

void
GrainBoundaryCavitation::initQpStatefulProperties()
{
  _a[_qp] = _a0[_qp];
  _b[_qp] = _b0[_qp];
  _D[_qp] = effectiveDamage();
  _ec[_qp] = 0;
  _delta_V[_qp] = 0;
  _delta_N[_qp] = 0;
}

void
GrainBoundaryCavitation::computeQpProperties()
{
  computeDamage();
  computeInterfaceTractionAndDerivatives();
  computeDamageDrivingForces();
  computeTimeStepLimit();
}

void
GrainBoundaryCavitation::computeDamage()
{
  const Real an = _a_old[_qp];
  const Real bn = _b_old[_qp];

  // If already fully damaged, don't bother updating
  if (MooseUtils::absoluteFuzzyEqual(_D_old[_qp], 1))
  {
    _a[_qp] = an;
    _b[_qp] = bn;
    _D[_qp] = 1;
  }
  // Otherwise, implicitly update the damage
  else
  {
    const Real h = cavityShapeFactor();

    // Update a
    if (_delta_V_old[_qp] > 0)
    {
      Real d = _delta_V_old[_qp] / 4 / M_PI / h;
      Real a2 = an * an;
      Real a3 = a2 * an;
      Real c =
          std::cbrt(3 * std::sqrt(12 * a3 * d + 81 * d * d) + 2 * a3 + 27 * d) / std::cbrt(2.0);
      _a[_qp] = (c + a2 / c + an) / 3;
    }
    else
      _a[_qp] = _a_old[_qp];

    // Update b
    if (_delta_N_old[_qp] > 1e-12)
    {
      Real m = M_PI / 2 * _delta_N_old[_qp];
      Real m2 = m * m;
      Real m3 = m2 * m;
      Real m4 = m3 * m;
      Real k = std::cbrt(std::sqrt(81 * bn * bn * m4 + 12 * m3) + 9 * bn * m2);
      _b[_qp] = k / std::cbrt(18.0) / m - std::cbrt(2.0 / 3.0) / k;
    }
    else
      _b[_qp] = _b_old[_qp];

    // Satisfy the constraints:
    // a <= b
    if (_a[_qp] > _b[_qp])
    {
      const Real da = _a[_qp] - an;
      const Real db = _b[_qp] - bn;
      const Real p = (bn - an) / (da - db);
      _a[_qp] = an + p * da;
      _b[_qp] = bn + p * db;
    }

    // Update damage
    _D[_qp] = effectiveDamage();
  }
}

void
GrainBoundaryCavitation::computeInterfaceTractionAndDerivatives()
{
  // Update normal traction
  const Real C_n = normalStiffness();
  const Real ju_n = normalSeparation();
  const auto [P, d_P_d_ju_n] = penetrationPenalty(ju_n);
  const Real T_n = P * C_n * ju_n;

  // Update tangential tractions
  const Real C_s = tangentialStiffness();
  const Real ju_s1 = _interface_displacement_jump[_qp](1);
  const Real ju_s2 = _interface_displacement_jump[_qp](2);
  const Real T_s1 = C_s * ju_s1;
  const Real T_s2 = C_s * ju_s2;

  _interface_traction[_qp] = RealVectorValue(T_n, T_s1, T_s2);
  _dinterface_traction_djump[_qp].zero();
  _dinterface_traction_djump[_qp](0, 0) = P * C_n + d_P_d_ju_n * C_n * ju_n;
  _dinterface_traction_djump[_qp](1, 1) = C_s;
  _dinterface_traction_djump[_qp](2, 2) = C_s;
}

void
GrainBoundaryCavitation::computeDamageDrivingForces()
{
  // Average stress
  const auto sigma = (_stress_primary[_qp] + _stress_secondary[_qp]) / 2;
  _sigma_h[_qp] = sigma.trace() / 3;
  _sigma_vm[_qp] = std::sqrt(3 * sigma.secondInvariant());
  _eta[_qp] =
      MooseUtils::absoluteFuzzyEqual(_sigma_vm[_qp], 0) ? 0 : _sigma_h[_qp] / _sigma_vm[_qp];

  // Average creep strain
  const auto Ec = (_creep_strain_primary[_qp] + _creep_strain_secondary[_qp]) / 2;
  _ec[_qp] = std::sqrt(2.0 / 3.0 * Ec.doubleContraction(Ec) + libMesh::TOLERANCE);
  _delta_ec[_qp] = _ec[_qp] - _ec_old[_qp];

  // Normal traction
  _Tn[_qp] = _interface_traction[_qp](0);

  // Number density (increment)
  const Real Nr = _FN[_qp] * std::pow(MathUtils::macaulay(_Tn[_qp]) / _T0[_qp], _gamma);
  _nucleating[_qp] = _nucleating_old[_qp];
  if (_nucleating_old[_qp] < 0.5)
    _nucleating[_qp] = MathUtils::heaviside(Nr * _ec[_qp] - _Nc[_qp]);
  _delta_N[_qp] = Nr * _delta_ec[_qp] * _nucleating_old[_qp];

  // Cavity volume increment
  _delta_V[_qp] = 0;
  if (_diffusion_growth)
  {
    const Real L =
        std::cbrt(_D_GB[_qp] * _sigma_vm[_qp] * _dt / (_delta_ec[_qp] + libMesh::TOLERANCE));
    const Real vf = MooseUtils::absoluteFuzzyEqual(L, 0)
                        ? std::pow(_a[_qp] / _b[_qp], 2.0)
                        : std::max(std::pow(_a[_qp] / _b[_qp], 2.0),
                                   std::pow(_a[_qp] / (_a[_qp] + 1.5 * L), 2.0));
    const Real q = -2.0 * std::log(vf) - (1.0 - vf) * (3.0 - vf);
    _delta_V[_qp] += 8.0 * M_PI * _D_GB[_qp] * _Tn[_qp] / q;
  }
  if (_creep_growth)
  {
    if (!_fixed_triaxiality_state)
    {
      const Real g =
          _sigma_h[_qp] > 0 ? std::log(3.0) - 2.0 / 3.0 : 2.0 * M_PI / 9.0 / std::sqrt(3.0);
      const Real alpha = 1.5 / _n;
      const Real beta = (_n - 1.0) * (_n + g) / _n / _n;
      const Real V = cavityVolume();
      _delta_V[_qp] += degradation() * 1.5 * _delta_ec[_qp] * V *
                       (std::abs(_eta[_qp]) > 1
                            ? MathUtils::sign(_eta[_qp]) * std::pow(alpha * _eta[_qp] + beta, _n)
                            : _eta[_qp] * std::pow(alpha + beta, _n));
    }
    else
    {
      if (_triaxiality_state == "LOW")
      {
        const Real g = 2.0 * M_PI / 9.0 / std::sqrt(3.0);
        const Real alpha = 1.5 / _n;
        const Real beta = (_n - 1.0) * (_n + g) / _n / _n;
        const Real V = cavityVolume();
        _delta_V[_qp] +=
            degradation() * 1.5 * _delta_ec[_qp] * V * (_eta[_qp] * std::pow(alpha + beta, _n));
      }
      else if (_triaxiality_state == "MEDIUM")
      {
        const Real g = std::log(3.0) - 2.0 / 3.0;
        const Real alpha = 1.5 / _n;
        const Real beta = (_n - 1.0) * (_n + g) / _n / _n;
        const Real V = cavityVolume();
        _delta_V[_qp] +=
            degradation() * 1.5 * _delta_ec[_qp] * V * (_eta[_qp] * std::pow(alpha + beta, _n));
      }
      else if (_triaxiality_state == "HIGH")
      {
        const Real g = std::log(3.0) - 2.0 / 3.0;
        const Real alpha = 1.5 / _n;
        const Real beta = (_n - 1.0) * (_n + g) / _n / _n;
        const Real V = cavityVolume();
        _delta_V[_qp] += degradation() * 1.5 * _delta_ec[_qp] * V *
                         (MathUtils::sign(_eta[_qp]) * std::pow(alpha * _eta[_qp] + beta, _n));
      }
    }
  }
}

Real
GrainBoundaryCavitation::cavityVolume() const
{
  return 4.0 / 3.0 * M_PI * std::pow(_a[_qp], 3.0) * cavityShapeFactor();
}

Real
GrainBoundaryCavitation::cavityShapeFactor() const
{
  return 1.0 / (1.0 + std::cos(_psi) - std::cos(_psi) / 2.0) / std::sin(_psi);
}

Real
GrainBoundaryCavitation::degradation() const
{
  return (1 - std::pow(_D[_qp], _p)) * (1 - _eps) + _eps;
}

Real
GrainBoundaryCavitation::effectiveDamage() const
{
  return _a[_qp] / _b[_qp];
}

Real
GrainBoundaryCavitation::normalSeparation() const
{
  // Macro normal separation
  Real ju_n_total = _interface_displacement_jump[_qp](0);

  // Micro normal separation
  Real V = cavityVolume();
  Real ju_n_micro = V / M_PI / _b[_qp] / _b[_qp];

  return ju_n_total - ju_n_micro;
}

std::pair<Real, Real>
GrainBoundaryCavitation::penetrationPenalty(Real ju_n) const
{
  const Real penalty = (1.0 + _P * std::pow(MathUtils::macaulay(-ju_n), 2.0));
  const Real d_penalty_d_ju_n = -2 * _P * MathUtils::macaulay(-ju_n);
  return {penalty, d_penalty_d_ju_n};
}

Real
GrainBoundaryCavitation::normalStiffness() const
{
  return degradation() * _E[_qp] / _w[_qp];
}

Real
GrainBoundaryCavitation::tangentialStiffness() const
{
  return degradation() * _G[_qp] / _w[_qp] / (1 - _G[_qp] / _w[_qp] / _eta_s[_qp] * _dt);
}

void
GrainBoundaryCavitation::computeTimeStepLimit()
{
  _dt_max[_qp] = _eta_s[_qp] / (_G[_qp] / _w[_qp]);
  Real delta_D = _D[_qp] - _D_old[_qp];
  if (delta_D > _delta_D_max)
    _dt_max[_qp] = std::min(_dt_max[_qp], _dt_cutback_factor * _dt / (delta_D / _delta_D_max));
}
