//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EffectiveStressTools.h"
#include "InterfaceEffectiveStress.h"
#include "InterfaceValueTools.h"

registerMooseObject("DeerApp", InterfaceEffectiveStress);

template <> InputParameters validParams<InterfaceEffectiveStress>() {
  InputParameters params = validParams<InterfaceMaterial>();
  params.addClassDescription("Compute an effective stress across an interface");
  params.addRequiredParam<MooseEnum>("effective_stress_type",
                                     EffectiveStressTools::scalarOptions(),
                                     "Type of effective stress to be computed");
  params.addParam<MooseEnum>("interface_value_type",
                             InterfaceValueTools::InterfaceAverageOptions(),
                             "Type of scalar output");
  params.addParam<std::vector<Real>>("params_vector",
                                     "vector of effective stress parameters");
  params.addParam<Point>("point1", Point(0, 0, 0),
                         "Start point for axis used to calculate some "
                         "cylindrical material tensor quantities");
  params.addRequiredParam<MaterialPropertyName>(
      "effective_stress_mp_name", "the name of the calcualte effective stress");
  params.addParam<Point>(
      "point2", Point(0, 1, 0),
      "End point for axis used to calculate some material tensor quantities");
  params.addParam<Point>("direction", Point(0, 0, 1), "Direction vector");
  return params;
}

InterfaceEffectiveStress::InterfaceEffectiveStress(
    const InputParameters &parameters)
    : InterfaceMaterial(parameters),
      _tensor_master(getMaterialPropertyByName<RankTwoTensor>("stress")),
      _tensor_slave(getNeighborMaterialPropertyByName<RankTwoTensor>("stress")),
      _effective_stress(declareProperty<Real>(
          getParam<MaterialPropertyName>("effective_stress_mp_name"))),
      _effective_stress_type(getParam<MooseEnum>("effective_stress_type")),
      _interface_value_type(parameters.get<MooseEnum>("interface_value_type")),
      _params_vector(getParam<std::vector<Real>>("params_vector")),
      _point1(parameters.get<Point>("point1")),
      _point2(parameters.get<Point>("point2")),
      _input_direction(parameters.get<Point>("direction") /
                       parameters.get<Point>("direction").norm()) {
  if (_effective_stress_type == 2 && _params_vector.size() != 1)
    mooseError("the huddleston effective stress requires the parameters b to "
               "be defined");
  if (_effective_stress_type == 3) {
    if (_params_vector.size() != 2)
      mooseError("the hayhurs effective stress requires the parameters alpha, "
                 "beta and gamma to be defined");
    _params_vector.push_back(1. - _params_vector[0] - _params_vector[1]);
  }
  if (_effective_stress_type == 6 && _params_vector.size() != 1)
    mooseError(
        "the RCCMRXMises effective stress requires the parameters alpha to "
        "be defined");
  if (_effective_stress_type == 7 && _params_vector.size() != 1)
    mooseError(
        "the RCCMRXTresca effective stress requires the parameters alpha to "
        "be defined");
}

void InterfaceEffectiveStress::computeQpProperties() {

  Real eff_stress_master = EffectiveStressTools::getQuantity(
      _tensor_master[_qp], _effective_stress_type, _params_vector, _point1,
      _point2, _input_direction);

  Real eff_stress_slave = EffectiveStressTools::getQuantity(
      _tensor_slave[_qp], _effective_stress_type, _params_vector, _point1,
      _point2, _input_direction);

  _effective_stress[_qp] = InterfaceValueTools::getQuantity(
      _interface_value_type, eff_stress_master, eff_stress_slave);
}

void InterfaceEffectiveStress::initQpStatefulProperties() {
  _effective_stress[_qp] = 0.0;
}
