#include "TimeDerivAuxKernel.h"

registerMooseObject("DeerApp", TimeDerivAuxKernel);

InputParameters
TimeDerivAuxKernel::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("coupled", "Coupled variable");

  return params;
}

TimeDerivAuxKernel::TimeDerivAuxKernel(const InputParameters & parameters)
  : AuxKernel(parameters),
    _coupled_new(coupledValue("coupled")),
    _coupled_old(coupledValueOld("coupled"))
{
}

Real
TimeDerivAuxKernel::computeValue()
{
  return (_coupled_new[_qp] - _coupled_old[_qp]) / _dt;
}
