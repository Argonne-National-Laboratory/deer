#include "AdiabaticHeating.h"

registerMooseObject("DeerApp", AdiabaticHeating);

InputParameters
AdiabaticHeating::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addClassDescription("Adds a heat source based on the value of a coupled variable.");

  params.addRequiredCoupledVar("heat", "Variable with the heating power");
  params.addParam<Real>("fraction", 1.0, "Optional value to multiply the input heat before adding into the system");

  return params;
}

AdiabaticHeating::AdiabaticHeating(const InputParameters & parameters) :
    Kernel(parameters),
    _heat(coupledValue("heat")),
    _fraction(getParam<Real>("fraction"))
{
}

Real
AdiabaticHeating::computeQpResidual()
{
  return _test[_i][_qp] * _fraction * _heat[_qp];
}
