//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BiLinearMixedModeTractionNewer.h"
#include "libmesh/utility.h"

registerMooseObject("DeerApp", BiLinearMixedModeTractionNewer);

InputParameters
BiLinearMixedModeTractionNewer::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addRequiredParam<Real>("penalty_stiffness", "Penalty stiffness.");
  params.addRequiredParam<Real>("GI_C", "Critical energy release rate in normal direction.");
  params.addRequiredParam<Real>("GII_C", "Critical energy release rate in shear direction.");
  params.addRequiredParam<Real>("normal_strength", "Tensile strength in normal direction.");
  params.addRequiredParam<Real>("shear_strength", "Tensile strength in shear direction.");
  params.addRequiredParam<Real>("eta", "The power law parameter.");
  MooseEnum criterion("POWER_LAW BK", "BK");
  params.addParam<MooseEnum>(
      "mixed_mode_criterion", criterion, "Option for mixed mode propagation criterion.");
  params.addClassDescription("Mixed mode bilinear traction separation law.");
  return params;
}

BiLinearMixedModeTractionNewer::BiLinearMixedModeTractionNewer(const InputParameters & parameters)
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
BiLinearMixedModeTractionNewer::computeInterfaceTractionAndDerivatives()
{
  // Calculate the mode-mixity -- lagged!
  _beta[_qp] = modeMixity(_interface_displacement_jump_old[_qp]);

  // Get the effective and threshold deltas
  Real delta_e = delta_eff(_interface_displacement_jump[_qp]);
  Real delta_i, delta_f;
  std::tie(delta_i, delta_f) = delta_threshold(_interface_displacement_jump[_qp]);

  // Calculate the prospective and actual damage
  Real d = damage(delta_e, delta_i, delta_f);
  if (d > _d_old[_qp])
    _d[_qp] = d;
  else 
    _d[_qp] = _d_old[_qp];

  // Calculate the stiffness
  _dinterface_traction_djump[_qp] = (1 - _d[_qp]) * _K * RankTwoTensor::Identity();
  if (_interface_displacement_jump[_qp](0) < 0)
    _dinterface_traction_djump[_qp](0,0) += _d[_qp] * _K;

  // Calculate the tractions
  _interface_traction[_qp] = _dinterface_traction_djump[_qp] * _interface_displacement_jump[_qp];

  // Account for the effect of damage on the derivative
  RealVectorValue d1 = -_K * _interface_displacement_jump[_qp];
  if (_interface_displacement_jump[_qp](0) < 0)
    d1(0) += _K * _interface_displacement_jump[_qp](0);

  RealVectorValue d2 = _interface_displacement_jump[_qp] / delta_e;
  if (_interface_displacement_jump[_qp](0) >= 0)
    d1(0) = 0;

  d2 *= (delta_f * delta_i) / (Utility::pow<2>(delta_e)*(delta_f - delta_i));

  RankTwoTensor extra;
  extra.vectorOuterProduct(d1, d2);

  _dinterface_traction_djump[_qp] += extra;
}

Real
BiLinearMixedModeTractionNewer::modeMixity(const RealVectorValue & delta) const
{
  if (delta(0) < libMesh::TOLERANCE * libMesh::TOLERANCE)
    return 0;
  const Real delta_shear = std::sqrt(delta(1) * delta(1) + delta(2) * delta(2));
  return delta_shear / delta(0);
}

void
BiLinearMixedModeTractionNewer::initQpStatefulProperties()
{
  CZMComputeLocalTractionTotalBase::initQpStatefulProperties();

  _d[_qp] = 0.0;
}

Real
BiLinearMixedModeTractionNewer::damage(const Real & delta,
                                     const Real & delta_init,
                                     const Real & delta_final) const
{
  if (delta < delta_init)
    return 0;
  if (delta > delta_final)
    return 1;

  return delta_final * (delta - delta_init) / delta / (delta_final - delta_init);
}

Real
BiLinearMixedModeTractionNewer::delta_eff(const RealVectorValue & delta) const
{
  Real delta_eff = Utility::pow<2>(delta(1)) + Utility::pow<2>(delta(2));
  if (delta(0) > 0.0)
    delta_eff += Utility::pow<2>(delta(0));
  return std::sqrt(delta_eff);
}

std::tuple<Real,Real>
BiLinearMixedModeTractionNewer::delta_threshold(const RealVectorValue & delta) const
{
  Real delta_o = _N / _K;
  Real delta_s = _S / _K;
  Real delta_init;

  if (delta(0) > 0.0)
    delta_init =  delta_o*delta_s*std::sqrt(1.0+Utility::pow<2>(_beta[_qp]) / 
                                     (Utility::pow<2>(delta_s) + 
                                      Utility::pow<2>(delta_o * _beta[_qp])));
  else
    delta_init = delta_s;

  Real delta_final;

  if (delta(0) > 0.0)
  {
    if (_criterion == MixedModeCriterion::BK)
      delta_final =  2.0/(_K*delta_init)*(_GI_C + (_GII_C - _GI_C))*std::pow(
              Utility::pow<2>(_beta[_qp])/(1+Utility::pow<2>(_beta[_qp])), _eta);
    else if (_criterion == MixedModeCriterion::POWER_LAW)
      delta_final = 2.0*(1.0+Utility::pow<2>(_beta[_qp]))/(_K*delta_init)*std::pow(
          std::pow(1.0/_GI_C,_eta)+std::pow(Utility::pow<2>(_beta[_qp])/_GII_C,_eta),-1.0/_eta);
  }
  else
    delta_final = std::sqrt(2) * 2 * _GII_C / _S;


  return std::make_tuple(delta_init, delta_final);

}
