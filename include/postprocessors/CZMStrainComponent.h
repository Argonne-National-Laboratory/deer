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
#include "CZMAreaRatioPostprocessor.h"

/**
 * This postprocessor computes the cohesive zone contribution to the total RVE
 * strain
 */

class CZMStrainComponent : public CZMAreaRatioPostprocessor
{
public:
  static InputParameters validParams();

  CZMStrainComponent(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  virtual Real computeStrainIntegral();

  /// the czm volumetric strain tensor
  const MaterialProperty<RankTwoTensor> & _tensor;

  /// the index i of the strain tensor
  const unsigned int _i;

  /// the index j of the strain tensor
  const unsigned int _j;

  /// the postprocessor value representing the undeformed RVE volume
  const PostprocessorValue & _initial_bulk_volume_pp;

  /// the czm strain value
  Real _normalized_strain_component;
};
