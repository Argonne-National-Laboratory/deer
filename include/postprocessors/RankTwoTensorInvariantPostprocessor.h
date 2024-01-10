//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"
#include "RankTwoScalarTools.h"

/**
 * Computes an invariant of a symmetric RankTwoTensor given a set of properly
 * named postprocessors that toghether represent a tensor.
 */
class RankTwoTensorInvariantPostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  RankTwoTensorInvariantPostprocessor(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override{};
  virtual void execute() override;
  /// we don't need to override finalize because we are working with
  /// postprocessors values. Therefore the values this PP will get, have already
  /// be summed/averaged among processes/threads.
  virtual PostprocessorValue getValue() const;

protected:
  const PostprocessorName _rank_two_tensor_base_name;
  /**
   * Determines the tensor invariant to be computed  using the
   * RankTwoScalarTools namespace, e.g., vonMisesStress L2norm, MaxPrincipal
   * eigenvalue, etc.
   */
  RankTwoScalarTools::InvariantType _invariant_type;

  /// The postprocessor values that toghether represent a rank two tensor
  std::vector<std::vector<const PostprocessorValue *>> _pps_values;

  /// the calculated invariant value to be returned
  PostprocessorValue _invariant;

  /// the map between postprocessor names and tensorial components
  const std::map<std::pair<int, int>, std::string> tensor_map = {{std::make_pair(0, 0), "xx"},
                                                                 {std::make_pair(1, 1), "yy"},
                                                                 {std::make_pair(2, 2), "zz"},
                                                                 {std::make_pair(0, 1), "xy"},
                                                                 {std::make_pair(0, 2), "xz"},
                                                                 {std::make_pair(1, 2), "yz"}};
};
