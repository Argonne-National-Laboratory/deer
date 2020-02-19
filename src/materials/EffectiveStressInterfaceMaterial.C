//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EffectiveStressInterfaceMaterial.h"
#include "EffectiveStressTools.h"
#include "InterfaceValueTools.h"

registerMooseObject("DeerApp", EffectiveStressInterfaceMaterial);

InputParameters EffectiveStressInterfaceMaterial::validParams() {
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription("Compute an effective stress across an interface");
  params.addRequiredParam<MooseEnum>("effective_stress_type",
                                     EffectiveStressTools::scalarOptions(),
                                     "Type of effective stress to be computed");
  params.addParam<MooseEnum>("interface_value_type",
                             InterfaceValueTools::InterfaceAverageOptions(),
                             "Type of scalar output");
  params.addParam<std::vector<Real>>("params_vector",
                                     "Vector of effective stress parameters");
  params.addRequiredParam<MaterialPropertyName>(
      "effective_stress_mp_name", "the name of the calcualte effective stress");
  params.addParam<bool>("stateful", false,
                        "If true make the material property stateful");
  return params;
}

EffectiveStressInterfaceMaterial::EffectiveStressInterfaceMaterial(
    const InputParameters &parameters)
    : InterfaceMaterial(parameters),
      _stress_master(getMaterialPropertyByName<RankTwoTensor>("stress")),
      _stress_slave(getNeighborMaterialPropertyByName<RankTwoTensor>("stress")),
      _effective_stress(declareProperty<Real>(
          getParam<MaterialPropertyName>("effective_stress_mp_name"))),
      _effective_stress_type(getParam<MooseEnum>("effective_stress_type")),
      _interface_value_type(parameters.get<MooseEnum>("interface_value_type")),
      _params_vector(getParam<std::vector<Real>>("params_vector")),
      _stateful(getParam<bool>("stateful")) {
  if (_effective_stress_type == 2 && _params_vector.size() != 1)
    mooseError("The Huddleston effective stress requires the parameters b to "
               "be defined");
  if (_effective_stress_type == 3) {
    if (_params_vector.size() != 2)
      mooseError("The Hayhurst effective stress requires the parameters alpha, "
                 "beta to be defined. Gamme is compute as 1-alpha-beta");
    _params_vector.push_back(1. - _params_vector[0] - _params_vector[1]);
  }
  if (_effective_stress_type == 6 && _params_vector.size() != 1)
    mooseError(
        "The RCCMRx-Mises effective stress requires the parameters alpha to "
        "be defined");
  if (_effective_stress_type == 7 && _params_vector.size() != 1)
    mooseError(
        "The RCCMRx-Tresca effective stress requires the parameters alpha to "
        "be defined");
}

void EffectiveStressInterfaceMaterial::computeQpProperties() {

  Real eff_stress_master = EffectiveStressTools::getQuantity(
      _stress_master[_qp], _effective_stress_type, _params_vector);

  Real eff_stress_slave = EffectiveStressTools::getQuantity(
      _stress_slave[_qp], _effective_stress_type, _params_vector);

  _effective_stress[_qp] = InterfaceValueTools::getQuantity(
      _interface_value_type, eff_stress_master, eff_stress_slave);
}

void EffectiveStressInterfaceMaterial::initQpStatefulProperties() {
  if (_stateful)
    _effective_stress[_qp] = 0.0;
}
