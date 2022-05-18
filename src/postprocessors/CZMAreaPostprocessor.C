//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMAreaPostprocessor.h"
#include "CohesiveZoneModelTools.h"

registerMooseObject("DeerApp", CZMAreaPostprocessor);

InputParameters
CZMAreaPostprocessor::validParams()
{
  InputParameters params = InterfaceIntegralPostprocessor::validParams();
  MooseEnum strainType("SMALL FINITE", "SMALL");
  params.addParam<MooseEnum>("strain", strainType, "Strain formulation");
  params.addParam<std::string>("base_name", "Material property base name");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addClassDescription("Computes the area of a cohesive zone based on "
                             "the provided strain formulation.");

  return params;
}

CZMAreaPostprocessor::CZMAreaPostprocessor(const InputParameters & parameters)
  : InterfaceIntegralPostprocessor(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _strain(getParam<MooseEnum>("strain").getEnum<Strain>()),
    _F_czm(_strain == Strain::Finite
               ? &getMaterialPropertyByName<RankTwoTensor>(_base_name + "F_czm")
               : nullptr)
{
}

Real
CZMAreaPostprocessor::computeQpIntegral()
{

  Real dadA = 1;

  if (_strain == Strain::Finite)
    dadA = CohesiveZoneModelTools::computeAreaRatio(
        (*_F_czm)[_qp].inverse().transpose(), (*_F_czm)[_qp].det(), _normals[_qp]);
  return dadA;
}
