//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "MooseEnum.h"
#include "MooseError.h"
#include "MooseTypes.h"
#include "RankTwoScalarTools.h"
#include "RankTwoTensor.h"
#include "libmesh/point.h"

namespace EffectiveStressTools {
/*
 * Return the scalar_type MooseEnum
 */
MooseEnum scalarOptions();

/*
 * Return an effective stress on the user specified
 * scalar_type
 * @param point1 The starting point of the rotation axis for a cylinderical
 * system
 * @param point2 The end point of the rotation axis
 * @param curr_point The point corresponding to the stress (pass in &
 * _q_point[_qp])
 * @param direction The direction vector in which the scalar stress value is
 * calculated point1 and point2 are required only for the cases of axialStress,
 * hoopStress and radialStress curr_point is required only for the cases of
 * hoopStress and radialStress direction is required only for
 * directionValueTensor for all other cases, these parameters will take the
 * default values
 */

template <typename T>
T huddleston(const RankTwoTensorTempl<T> &stress, const Real &b) {

  Real I1 = stress.trace();
  Real svm = RankTwoScalarTools::vonMisesStress(stress);

  if (svm == 0)
    return 0;

  Real SS = std::sqrt(svm * svm + stress.secondInvariant());

  if (SS == 0)
    return 0;

  return svm * std::exp(b * (I1 / SS));
}

template <typename T>
T Hayhurst(const RankTwoTensorTempl<T> &stress,
           const std::vector<Real> &params_vector, Point &direction) {

  Real I1 = stress.trace();
  Real svm = RankTwoScalarTools::vonMisesStress(stress);
  Real S1_positive =
      std::max(RankTwoScalarTools::maxPrincipal(stress, direction), 0.0);
  return params_vector[0] * S1_positive + params_vector[1] * I1 +
         params_vector[2] * svm;
}

template <typename T>
T RCCMRXMises(const RankTwoTensorTempl<T> &stress, const Real &alpha) {
  Real svm = RankTwoScalarTools::vonMisesStress(stress);
  return alpha * stress.trace() + (1. - alpha) * svm;
}

template <typename T>
T RCCMRXTresca(const RankTwoTensorTempl<T> &stress, const Real &alpha) {

  Real tresca = RankTwoScalarTools::stressIntensity(stress);
  return alpha * stress.trace() + (1. - alpha) * tresca;
}

template <typename T>
T getQuantity(const RankTwoTensorTempl<T> &stress, const MooseEnum &scalar_type,
              const std::vector<Real> &params_vector) {
  Point direction(1, 0, 0);
  switch (scalar_type) {
  case 0:
    return RankTwoScalarTools::vonMisesStress(stress);
  case 1:
    return RankTwoScalarTools::hydrostatic(stress);
  case 2:
    return huddleston(stress, params_vector[0]);
  case 3:
    return Hayhurst(stress, params_vector, direction);
  case 4:
    return RankTwoScalarTools::maxPrincipal(stress, direction);
  case 5:
    return RankTwoScalarTools::stressIntensity(stress);
  case 6:
    return RCCMRXMises(stress, params_vector[0]);
  case 7:
    return RCCMRXTresca(stress, params_vector[0]);
  case 8:
    return std::max(RankTwoScalarTools::maxPrincipal(stress, direction),
                    RankTwoScalarTools::vonMisesStress(stress));
  default:
    mooseError("RankTwoScalarAux Error: Pass valid scalar type - " +
               scalarOptions().getRawNames());
  }
}

} // namespace EffectiveStressTools
