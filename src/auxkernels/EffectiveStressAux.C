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

defineLegacyParams(EffectiveStressAux);

InputParameters EffectiveStressAux::validParams() {
  InputParameters params = NodalPatchRecovery::validParams();
  params.addClassDescription("Compute an effective stress");
  params.addParam<MooseEnum>("effective_stress_type",
                             EffectiveStressTools::scalarOptions(),
                             "Type of scalar output");
  params.addParam<std::vector<Real>>("params_vector",
                                     "vector of effective stress parameters");
  params.addParam<Point>("point1", Point(0, 0, 0),
                         "Start point for axis used to calculate some "
                         "cylindrical material tensor quantities");
  params.addParam<Point>(
      "point2", Point(0, 1, 0),
      "End point for axis used to calculate some material tensor quantities");
  params.addParam<Point>("direction", Point(0, 0, 1), "Direction vector");
  return params;
}

EffectiveStressAux::EffectiveStressAux(const InputParameters &parameters)
    : NodalPatchRecovery(parameters),
      _tensor(getMaterialProperty<RankTwoTensor>("stress")),
      _scalar_type(getParam<MooseEnum>("effective_stress_type")),
      _params_vector(getParam<std::vector<Real>>("params_vector")),
      _point1(parameters.get<Point>("point1")),
      _point2(parameters.get<Point>("point2")),
      _input_direction(parameters.get<Point>("direction") /
                       parameters.get<Point>("direction").norm()) {

  if (_scalar_type == 2 && _params_vector.size() != 1)
    mooseError("the huddleston effective stress requires the parameters b to "
               "be defined");
  if (_scalar_type == 3) {
    if (_params_vector.size() != 2)
      mooseError("the hayhurs effective stress requires the parameters alpha, "
                 "beta and gamma to be defined");
    _params_vector.push_back(1. - _params_vector[0] - _params_vector[1]);
  }
  if (_scalar_type == 6 && _params_vector.size() != 1)
    mooseError(
        "the RCCMRXMises effective stress requires the parameters alpha to "
        "be defined");
  if (_scalar_type == 7 && _params_vector.size() != 1)
    mooseError(
        "the RCCMRXTresca effective stress requires the parameters alpha to "
        "be defined");
}

Real EffectiveStressAux::computeValue() {

  return EffectiveStressTools::getQuantity(_tensor[_qp], _scalar_type,
                                           _params_vector, _point1, _point2,
                                           _input_direction);
}
