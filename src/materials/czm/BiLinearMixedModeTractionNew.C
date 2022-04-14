//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BiLinearMixedModeTractionNew.h"
#include "libmesh/utility.h"

registerMooseObject("DeerApp", BiLinearMixedModeTractionNew);

InputParameters
BiLinearMixedModeTractionNew::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addRequiredParam<Real>("penalty_stiffness", "Penalty stiffness.");
  params.addRequiredParam<Real>("GI_C", "Critical energy release rate in normal direction.");
  params.addRequiredParam<Real>("GII_C", "Critical energy release rate in shear direction.");
  params.addRequiredParam<Real>("normal_strength", "Tensile strength in normal direction.");
  params.addRequiredParam<Real>("shear_strength", "Tensile strength in shear direction.");
  params.addRequiredParam<Real>("eta", "The power law parameter.");
  params.addParam<bool>(
      "lag_seperation_state", false, "Lag seperation sate: use their old values.");
  MooseEnum criterion("POWER_LAW BK", "BK");
  params.addParam<MooseEnum>(
      "mixed_mode_criterion", criterion, "Option for mixed mode propagation criterion.");
  params.addClassDescription("Mixed mode bilinear traction separation law.");
  return params;
}

BiLinearMixedModeTractionNew::BiLinearMixedModeTractionNew(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _K(getParam<Real>("penalty_stiffness")),
    _d(declareProperty<Real>("damage")),
    _d_old(getMaterialPropertyOld<Real>("damage")),
    _interface_displacement_jump_old(
        getMaterialPropertyOld<RealVectorValue>("interface_displacement_jump")),
    _GI_C(getParam<Real>("GI_C")),
    _GII_C(getParam<Real>("GII_C")),
    _N(getParam<Real>("normal_strength")),
    _S(getParam<Real>("shear_strength")),
    _eta(getParam<Real>("eta")),
    _beta(declareProperty<Real>("mode_mixity_ratio")),
    _criterion(getParam<MooseEnum>("mixed_mode_criterion").getEnum<MixedModeCriterion>())
{
}

void
BiLinearMixedModeTractionNew::computeInterfaceTractionAndDerivatives()
{
  _interface_traction[_qp] = computeTraction();
  _dinterface_traction_djump[_qp] = computeTractionDerivatives();
}

void
BiLinearMixedModeTractionNew::initQpStatefulProperties()
{
  CZMComputeLocalTractionTotalBase::initQpStatefulProperties();

  _d[_qp] = 0.0;
}

Real
BiLinearMixedModeTractionNew::modeMixity(const RealVectorValue & delta)
{
  if (delta(0) < libMesh::TOLERANCE * libMesh::TOLERANCE)
    return 0;
  const Real delta_shear = std::sqrt(delta(1) * delta(1) + delta(2) * delta(2));
  return delta_shear / delta(0);
}

Real
BiLinearMixedModeTractionNew::damage(const Real & delta,
                                     const Real & delta_init,
                                     const Real & delta_final)
{
  if (delta < delta_init)
    return 0;
  if (delta > delta_final)
    return 1;

  return delta_final * (delta - delta_init) / delta / (delta_final - delta_init);
}

RealVectorValue
BiLinearMixedModeTractionNew::computeTraction()
{
  const Real delta_normal0 = _N / _K;
  const Real delta_shear0 = _S / _K;

  const RealVectorValue delta = _interface_displacement_jump[_qp];
  const RealVectorValue delta_old = _interface_displacement_jump_old[_qp];

  // Compute mode mixity
  const Real beta = modeMixity(delta_old);

  // Compute relative displacement jump at damage initiation
  Real delta_init = 0;
  if (delta(0) > 0)
    delta_init = delta_normal0 * delta_shear0 *
                 std::sqrt((1 + beta * beta) /
                           (delta_shear0 * delta_shear0 + Utility::pow<2>(beta * delta_normal0)));
  else
    delta_init = delta_shear0;

  // Compute relative displacement jump at full degradation
  Real delta_final = 0;
  if (delta(0) > 0)
  {
    if (_criterion == MixedModeCriterion::BK)
      delta_final = 2 / _K / delta_init *
                    (_GI_C + (_GII_C - _GI_C) * std::pow(beta * beta / (1 + beta * beta), _eta));
    else if (_criterion == MixedModeCriterion::POWER_LAW)
      delta_final =
          (2 + 2 * beta * beta) / _K / delta_init *
          std::pow(std::pow(1 / _GI_C, _eta) + std::pow(beta * beta / _GII_C, _eta), -1.0 / _eta);
  }
  else
    delta_final = std::sqrt(2) * 2 * _GII_C / _S;

  // Compute mixed mode relative displacement
  const Real delta_m = std::sqrt(Utility::pow<2>(delta_old(1)) + Utility::pow<2>(delta_old(2)) +
                                 Utility::pow<2>(std::max(delta_old(0), 0.)));

  // Compute damage variable
  const Real d = damage(delta_m, delta_init, delta_final);

  // Irreversibility
  if (d > _d_old[_qp])
    _d[_qp] = d;
  else
    _d[_qp] = _d_old[_qp];

  // Split displacement jump into active and inactive parts
  const RealVectorValue delta_active(std::max(delta(0), 0.), delta(1), delta(2));
  const RealVectorValue delta_inactive(std::min(delta(0), 0.), 0, 0);

  return (1 - _d[_qp]) * _K * delta_active + _K * delta_inactive;
}

RankTwoTensor
BiLinearMixedModeTractionNew::computeTractionDerivatives()
{
  const RealVectorValue delta = _interface_displacement_jump[_qp];
  RankTwoTensor ddelta_active_ddelta, ddelta_inactive_ddelta;
  ddelta_active_ddelta.fillFromInputVector({delta(0) > 0 ? 1. : 0., 1., 1.});
  ddelta_inactive_ddelta.fillFromInputVector({delta(0) < 0 ? 1. : 0., 0., 0.});
  return (1 - _d[_qp]) * _K * ddelta_active_ddelta + _K * ddelta_inactive_ddelta;
}
