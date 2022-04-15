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
  params.addParam<bool>("lag_mode_mixity", false, "Lag the mode mixity by one "
                        "step, can help with convergence");
  params.addParam<bool>("lag_damage", false, "Lag the damage by one step, "
                        "can help with convergence");
  params.addClassDescription("Mixed mode bilinear traction separation law.");
  return params;
}

BiLinearMixedModeTractionNewer::BiLinearMixedModeTractionNewer(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _interface_displacement_jump_old(
        getMaterialPropertyOld<RealVectorValue>("interface_displacement_jump")),
    _K(getParam<Real>("penalty_stiffness")),
    _d(declareProperty<Real>("damage")),
    _d_old(getMaterialPropertyOld<Real>("damage")),
    _GI_C(getParam<Real>("GI_C")),
    _GII_C(getParam<Real>("GII_C")),
    _N(getParam<Real>("normal_strength")),
    _S(getParam<Real>("shear_strength")),
    _eta(getParam<Real>("eta")),
    _criterion(getParam<MooseEnum>("mixed_mode_criterion").getEnum<MixedModeCriterion>()),
    _lag_beta(getParam<bool>("lag_mode_mixity")),
    _lag_damage(getParam<bool>("lag_damage"))
{
}

void
BiLinearMixedModeTractionNewer::computeInterfaceTractionAndDerivatives()
{
  std::tie(_interface_traction[_qp], 
           _dinterface_traction_djump[_qp]) = 
      updateState(_interface_displacement_jump[_qp]);
}

void
BiLinearMixedModeTractionNewer::initQpStatefulProperties()
{
  CZMComputeLocalTractionTotalBase::initQpStatefulProperties();

  _d[_qp] = 0.0;
}

std::tuple<RealVectorValue, RankTwoTensor>
BiLinearMixedModeTractionNewer::updateState(const RealVectorValue & delta)
{
  RealVectorValue traction;
  RankTwoTensor jacobian;

  // Calculate the mode mixity and derivative
  Real beta;
  RealVectorValue dbeta;
  if (_lag_beta)
    std::tie(beta, dbeta) = modeMixity(_interface_displacement_jump_old[_qp]);
  else
    std::tie(beta, dbeta) = modeMixity(delta);
  
  // Calculate the effective jump and derivative
  auto [de, dde] = deltaEffective(delta);

  // Calculate the initial and final jumps and derivatives
  auto [di, ddi, df, ddf] = deltaThreshold(delta, beta);

  // Calculate the potential damage and associated derivatives
  auto [dm, dme, dmi, dmf] = damage(de, di, df);

  // Determine if we're going to accept the damage
  bool accept = dm > _d_old[_qp];
  if (accept)
    _d[_qp] = dm;
  else
    _d[_qp] = _d_old[_qp];
  
  Real duse = _d[_qp];
  if (_lag_damage) {
    duse = _d_old[_qp];
    accept = false;
  }
  
  // Calculate the operator
  bool comp = delta(0) < 0;
  jacobian = RankTwoTensor::Identity() * (1 - duse) * _K;
  if (comp)
    jacobian(0,0) += duse * _K;

  // Calculate the updated traction
  traction = jacobian * delta;

  // If we accepted the step then get the extra "d_damage"
  // terms for the derivative
  if (accept) {
    RealVectorValue v1 = -_K*delta;
    if (comp)
      v1(0) += _K * delta(0);

    RealVectorValue v2 = dme  * dde;
    if (!_lag_beta)
      v2 += (dmi * ddi + dmf * ddf) * dbeta;
    
    // This feels clunky
    RankTwoTensor extra;
    extra.vectorOuterProduct(v1, v2);

    // Add the extra term
    jacobian += extra;
  }

  return {traction, jacobian};
}

std::tuple<Real, RealVectorValue>
BiLinearMixedModeTractionNewer::modeMixity(const RealVectorValue & delta) const
{
  Real beta;
  RealVectorValue dbeta;

  if (delta(0) >= libMesh::TOLERANCE * libMesh::TOLERANCE) {
    beta = std::sqrt(delta(1) * delta(1) + delta(2) * delta(2)) / delta(0);
    dbeta(0) = -beta / delta(0);
    dbeta(1) = delta(1) / (beta * delta(0) * delta(0));
    dbeta(2) = delta(2) / (beta * delta(0) * delta(0));
  }

  return {beta, dbeta};
}

std::tuple<Real,Real,Real,Real>
BiLinearMixedModeTractionNewer::damage(const Real & delta,
                                     const Real & delta_init,
                                     const Real & delta_final) const
{
  Real d, d_e, d_o, d_f;

  if (delta < delta_init) {
    d = 0.0;
    d_e = 0.0;
    d_o = 0.0;
    d_f = 0.0;
  }
  else if (delta <= delta_final) {
    d = delta_final * (delta - delta_init) / (delta * (delta_final -
                                                       delta_init));
    d_e = delta_final * delta_init / (delta * delta * (delta_final -
                                                       delta_init));
    d_o = delta_final * (delta - delta_final) / (delta *
                                                 Utility::pow<2>(delta_final -
                                                                 delta_init));
    d_f = delta_init * (delta_init - delta) / (delta *
                                               Utility::pow<2>(delta_final -
                                                               delta_init));
  }
  else {
    d = 1.0;
    d_e = 0.0;
    d_o = 0.0;
    d_f = 0.0;
  }

  return {d,d_e,d_o,d_f};
}

std::tuple<Real, RealVectorValue>
BiLinearMixedModeTractionNewer::deltaEffective(const RealVectorValue & delta) const
{
  Real delta_m = std::sqrt(Utility::pow<2>(delta(1)) + Utility::pow<2>(delta(2)) +
                         Utility::pow<2>(std::max(delta(0), 0.0)));
  RealVectorValue ddelta = delta / delta_m;
  if (delta(0) < 0.0) ddelta(0) = 0.0;
  return {delta_m, ddelta};
}

std::tuple<Real,Real,Real,Real>
BiLinearMixedModeTractionNewer::deltaThreshold(const RealVectorValue & delta,
                                                const Real & beta) const
{
  Real d1 = _N / _K;
  Real ds = _S / _K;

  Real d_o, dd_o, d_f, dd_f;
  if (delta(0) > 0) {
    // The onset jump and derivative
    d_o = d1 * ds * std::sqrt((1.0 + beta * beta) / (ds * ds + beta * d1 * beta
                                                     * d1));
    dd_o = beta * Utility::pow<2>(d1*ds) * 
        (ds*ds - d1*d1) / (d_o * 
                           Utility::pow<2>(beta*beta*d1*d1+ds*ds));

    // The final jump and derivative
    if (_criterion == MixedModeCriterion::BK) {
      d_f = 2.0 / (_K * d_o) * (_GI_C + (_GII_C -
                                        _GI_C)*std::pow(beta*beta/(1+beta*beta),
                                                       _eta));
      dd_f = 4.0 * std::pow(beta*beta/(1+beta*beta), _eta) * (
          _GII_C - _GI_C) * _eta / ((beta*beta*beta+beta) * _K * d_o) - 
          d_f / d_o * dd_o;
    }
    else if (_criterion == MixedModeCriterion::POWER_LAW) {
      d_f = 2 * (1 + beta*beta) / (_K * d_o) * 
          std::pow(std::pow(1.0/_GI_C, _eta) + 
                   std::pow(beta*beta / _GII_C, _eta), -1.0 / _eta);
      dd_f = 2.0 * d_f * (beta * beta * std::pow(1.0 / _GI_C, _eta) -
                          std::pow(beta*beta / _GII_C, _eta))*
          std::pow(std::pow(1.0/_GI_C, _eta) + 
                   std::pow(beta * beta / _GII_C, _eta), -1.0) / 
          (beta * (1.0 + beta*beta)) - d_f / d_o * dd_o;
    }
  }
  else {
    // The onset jump and derivative
    d_o = ds;
    dd_o = 0.0;

    // The final jump and derivative
    d_f = 2.0 * std::sqrt(2.0) * _GII_C / _S;
    dd_f = 0.0;
  }

  return {d_o,dd_o,d_f,dd_f};
}
