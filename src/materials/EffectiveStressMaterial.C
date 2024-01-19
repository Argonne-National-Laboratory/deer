//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EffectiveStressMaterial.h"
#include "EffectiveStressTools.h"
#include "InterfaceValueTools.h"

registerMooseObject("DeerApp", EffectiveStressMaterial);

InputParameters
EffectiveStressMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute an effective at each quadrature point");
  params.addRequiredParam<MooseEnum>("effective_stress_type",
                                     EffectiveStressTools::scalarOptions(),
                                     "Type of effective stress to be computed");
  params.addParam<std::vector<Real>>("params_vector", {}, "Vector of effective stress parameters");
  params.addRequiredParam<MaterialPropertyName>("effective_stress_mp_name",
                                                "The name of the new  material_property");
  params.addParam<bool>("stateful", false, "If true make the material property stateful");

  return params;
}

EffectiveStressMaterial::EffectiveStressMaterial(const InputParameters & parameters)
  : Material(parameters),
    _stress(getMaterialPropertyByName<RankTwoTensor>("cauchy_stress")),
    _effective_stress(
        declareProperty<Real>(getParam<MaterialPropertyName>("effective_stress_mp_name"))),
    _effective_stress_type(getParam<MooseEnum>("effective_stress_type")),
    _params_vector(getParam<std::vector<Real>>("params_vector")),
    _stateful(getParam<bool>("stateful"))
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

void
EffectiveStressMaterial::computeQpProperties()
{

  _effective_stress[_qp] =
      EffectiveStressTools::getQuantity(_stress[_qp], _effective_stress_type, _params_vector);
}

void
EffectiveStressMaterial::initQpStatefulProperties()
{
  if (_stateful)
    _effective_stress[_qp] = 0.0;
}
