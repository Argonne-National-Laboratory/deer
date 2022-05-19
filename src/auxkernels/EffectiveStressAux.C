//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EffectiveStressAux.h"
#include "EffectiveStressTools.h"

registerMooseObject("DeerApp", EffectiveStressAux);

InputParameters
EffectiveStressAux::validParams()
{
  InputParameters params = NodalPatchRecovery::validParams();
  params.addClassDescription("Compute an effective stress");
  params.addParam<MooseEnum>(
      "effective_stress_type", EffectiveStressTools::scalarOptions(), "Type of scalar output");
  params.addParam<std::vector<Real>>("params_vector", "Vector of effective stress parameters");
  return params;
}

EffectiveStressAux::EffectiveStressAux(const InputParameters & parameters)
  : NodalPatchRecovery(parameters),
    _tensor(getMaterialProperty<RankTwoTensor>("cauchy_stress")),
    _effective_stress_type(getParam<MooseEnum>("effective_stress_type")),
    _params_vector(getParam<std::vector<Real>>("params_vector"))
{
  if (_effective_stress_type == 2 && _params_vector.size() != 1)
    mooseError("The Huddleston effective stress requires the parameters b to "
               "be defined");
  if (_effective_stress_type == 3)
  {
    if (_params_vector.size() != 2)
      mooseError("The Hayhurst effective stress requires the parameters alpha, "
                 "beta to be defined. Gamme is compute as 1-alpha-beta");
    _params_vector.push_back(1. - _params_vector[0] - _params_vector[1]);
  }
  if (_effective_stress_type == 6 && _params_vector.size() != 1)
    mooseError("The RCCMRx-Mises effective stress requires the parameters alpha to "
               "be defined");
  if (_effective_stress_type == 7 && _params_vector.size() != 1)
    mooseError("The RCCMRx-Tresca effective stress requires the parameters alpha to "
               "be defined");
}

Real
EffectiveStressAux::computeValue()
{

  return EffectiveStressTools::getQuantity(_tensor[_qp], _effective_stress_type, _params_vector);
}
