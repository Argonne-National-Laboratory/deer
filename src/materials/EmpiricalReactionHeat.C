#include "EmpiricalReactionHeat.h"

registerMooseObject("DeerApp", EmpiricalReactionHeat);

InputParameters
EmpiricalReactionHeat::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<Real>("total_heat", "Total heat per volume produced by the reaction");
  params.addRequiredParam<Real>("start_temperature", "Temperature at which to start the reaction");
  params.addRequiredParam<Real>("time_constant", "Characteristic reaction time");

  params.addRequiredCoupledVar("temperature", "Coupled temperature");  

  return params;
}

EmpiricalReactionHeat::EmpiricalReactionHeat(const InputParameters & parameters)
  : Material(parameters),
    _W0(getParam<Real>("total_heat")),
    _T0(getParam<Real>("start_temperature")),
    _tau(getParam<Real>("time_constant")),
    _active(declareProperty<bool>("reaction_active")),
    _active_old(getMaterialPropertyOld<bool>("reaction_active")),
    _power(declareProperty<Real>("reaction_power")),
    _power_old(getMaterialPropertyOld<Real>("reaction_power")),
    _temperature(coupledValue("temperature"))
{

}

void
EmpiricalReactionHeat::computeQpProperties()
{
  if (!_active_old[_qp] && (_temperature[_qp] >= _T0)) {
    _active[_qp] = true;
    _power[_qp] = _W0 / _tau;
  }
  else if (_active_old[_qp]) {
    _active[_qp] = true;
    _power[_qp] = _power_old[_qp] * std::exp(-_dt/_tau);
  }
  else {
    _active[_qp] = false;
    _power[_qp] = 0.0;
  }
}

void
EmpiricalReactionHeat::initQpStatefulProperties()
{
  _active[_qp] = false;
  _power[_qp] = 0.0;
}
