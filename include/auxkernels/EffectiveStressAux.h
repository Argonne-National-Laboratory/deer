//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NodalPatchRecovery.h"
#include "RankTwoTensor.h"

class EffectiveStressAux;

template <> InputParameters validParams<EffectiveStressAux>();

/**
 * EffectiveStressAux uses the namespace EffectiveStressTools to compute scalar
 * values from Rank-2 tensors.
 */
class EffectiveStressAux : public NodalPatchRecovery {
public:
  static InputParameters validParams();

  EffectiveStressAux(const InputParameters &parameters);

protected:
  virtual Real computeValue();

  const MaterialProperty<RankTwoTensor> &_tensor;

  /**
   * Determines the information to be extracted from the tensor by using the
   * EffectiveStressTools namespace, e.g., vonMisesStressL2norm, MaxPrincipal
   * eigenvalue, etc.
   */
  MooseEnum _scalar_type;
  std::vector<Real> _params_vector;

  const Point _point1;
  const Point _point2;
  Point _input_direction;
};
