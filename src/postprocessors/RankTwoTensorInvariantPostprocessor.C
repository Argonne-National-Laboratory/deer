//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RankTwoTensorInvariantPostprocessor.h"

registerMooseObject("DeerApp", RankTwoTensorInvariantPostprocessor);

InputParameters
RankTwoTensorInvariantPostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription("Computes an invariant of a symmetric rank two "
                             "tensor given a set of postprocessors"
                             " representing tensorial components. This class "
                             "leverages MOOSE RankTwoScalarTools::invariantOptions()");
  params.addRequiredParam<PostprocessorName>(
      "rank_two_tensor_base_name",
      "The base name of the rank_two_tensor postprocessors from which the "
      "invariant shall to be computed.");
  MooseEnum mixedInvariants("VonMisesStress EffectiveStrain Hydrostatic L2norm "
                            "VolumetricStrain FirstInvariant "
                            "SecondInvariant ThirdInvariant TriaxialityStress "
                            "MaxShear StressIntensity MaxPrincipal "
                            "MidPrincipal MinPrincipal");

  params.addParam<MooseEnum>("invariant", mixedInvariants, "Type of invariant output");
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
  return params;
}

RankTwoTensorInvariantPostprocessor::RankTwoTensorInvariantPostprocessor(
    const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _rank_two_tensor_base_name(getParam<PostprocessorName>("rank_two_tensor_base_name")),
    _invariant_type(
        getParam<MooseEnum>("invariant").template getEnum<RankTwoScalarTools::InvariantType>())
{
}

void
RankTwoTensorInvariantPostprocessor::initialSetup()
{
  _pps_values.resize(3, std::vector<const PostprocessorValue *>(3));
  for (auto entry : tensor_map)
  {
    const int i = entry.first.first;
    const int j = entry.first.second;
    _pps_values[i][j] =
        (&getPostprocessorValueByName(_rank_two_tensor_base_name + "_" + entry.second));
    if (i != j)
      _pps_values[j][i] =
          (&getPostprocessorValueByName(_rank_two_tensor_base_name + "_" + entry.second));
  }
}

void
RankTwoTensorInvariantPostprocessor::execute()
{
  RankTwoTensor tensor;
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      tensor(i, j) = (*_pps_values[i][j]);
  switch (_invariant_type)
  {
    case RankTwoScalarTools::InvariantType::EffectiveStrain:
      _invariant = RankTwoScalarTools::effectiveStrain(tensor);
      break;
    default:
      _invariant = RankTwoScalarTools::getInvariant(tensor, _invariant_type);
      break;
  }
}

PostprocessorValue
RankTwoTensorInvariantPostprocessor::getValue() const
{
  return _invariant;
}
