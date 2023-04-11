//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ScalarPrincipalValues.h"

#include "RankTwoTensor.h"

registerMooseObject("DeerApp", ScalarPrincipalValues);

InputParameters
ScalarPrincipalValues::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredCoupledVar("scalar_variable", "Name of variable");
  params.addParam<size_t>("rank", 0, "Which principal component to return");
  return params;
}

ScalarPrincipalValues::ScalarPrincipalValues(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _scalar_var(coupledScalarValue("scalar_variable")),
    _scalar_order(coupledScalarOrder("scalar_variable")),
    _rank(getParam<size_t>("rank"))
{
  // Could add general case later
  if (_scalar_order != 6)
    mooseError("Coupled scalar variable needs to have order 6");

  if (_rank < 0)
    mooseError("Cannot request a rank less than 0");

  if (_rank > 2)
    mooseError("Cannot request a rank greater than 2");
}

void
ScalarPrincipalValues::initialize()
{
  _value = 0.0;
}

void
ScalarPrincipalValues::execute()
{
  RankTwoTensor val(_scalar_var[0],
                    _scalar_var[3],
                    _scalar_var[5],
                    _scalar_var[4],
                    _scalar_var[2],
                    _scalar_var[1]);

  std::vector<Real> vals;
  val.symmetricEigenvalues(vals);

  _value = vals[2 - _rank];
}

PostprocessorValue
ScalarPrincipalValues::getValue()
{
  return _value;
}
