#include "QuasiContactBC.h"

registerMooseObject("DeerApp", QuasiContactBC);

InputParameters
QuasiContactBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();

  params.addRequiredParam<FunctionName>("penalty_function", "The value of the penalty parameter");
  params.addParam<Real>("value", 0.0, "Value to enforce");

  return params;
}

QuasiContactBC::QuasiContactBC(const InputParameters & parameters)
  : IntegratedBC(parameters), _func(getFunction("penalty_function")), _v(getParam<Real>("value"))
{
}

Real
QuasiContactBC::computeQpResidual()
{
  return _func.value(_t, _q_point[_qp]) * _test[_i][_qp] * (-_v + _u[_qp]);
}

Real
QuasiContactBC::computeQpJacobian()
{
  return _func.value(_t, _q_point[_qp]) * _phi[_j][_qp] * _test[_i][_qp];
}
