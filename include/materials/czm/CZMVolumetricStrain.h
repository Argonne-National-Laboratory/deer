//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceMaterial.h"

/**
 * This material computes the interface strain tensor and separate it into its
 * normal and tangential contribution. If the large_kinematic flag is true
 * (default) area changes and large rotation are used to properly adjust the
 * resulting strains.
 */
class CZMVolumetricStrain : public InterfaceMaterial {

public:
  static InputParameters validParams();
  CZMVolumetricStrain(const InputParameters &parameters);

protected:
  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

  /// normal to the interface
  const MooseArray<Point> &_normals;

  /// Base name of the material system
  const std::string _base_name;

  /// the displacement jump in global and interface coordiantes
  ///@{
  const MaterialProperty<RealVectorValue> &_displacement_jump_global;
  ///@}

  /// the material property defining the czm normal
  const MaterialProperty<RankTwoTensor> &_czm_total_rotation;

  /// interface strains and strains rate
  ///@{
  MaterialProperty<RankTwoTensor> &_czm_total_strain;
  MaterialProperty<RankTwoTensor> &_czm_normal_strain;
  MaterialProperty<RankTwoTensor> &_czm_sliding_strain;
  ///@}

  /// ratio between deformed and undeformed area
  MaterialProperty<Real> &_dadA_mp;

  /// strain formulation
  enum class Strain { Small, Finite } _strain;

  const MaterialProperty<RankTwoTensor> *_F_czm;

  /// method computing the interface strain contribution
  void computeInterfaceStrain();
};
