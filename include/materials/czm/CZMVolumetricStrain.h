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
  virtual void initialSetup() override;

  virtual void computeQpProperties() override;

  /// normal to the interface
  const MooseArray<Point> &_normals;

  /// number of displacement components
  const unsigned int _ndisp;

  /// the coupled displacement and neighbor displacement values
  ///@{
  std::vector<const VariableValue *> _disp;
  std::vector<const VariableValue *> _disp_neighbor;
  std::vector<MooseVariable *> _disp_vars;
  std::vector<const VariableValue *> _disp_old;
  std::vector<const VariableValue *> _disp_neighbor_old;
  ///@}

  /// the coupled displacement and neighbor displacement gradient
  ///@{
  std::vector<const VariableGradient *> _grad_disp;
  std::vector<const VariableGradient *> _grad_disp_neighbor;
  std::vector<const VariableGradient *> _grad_disp_old;
  std::vector<const VariableGradient *> _grad_disp_neighbor_old;
  ///@}

  /// interface strains
  ///@{
  MaterialProperty<RankTwoTensor> &_czm_total_strain_rate;
  MaterialProperty<RankTwoTensor> &_czm_normal_strain_rate;
  MaterialProperty<RankTwoTensor> &_czm_sliding_strain_rate;
  ///@}

  /// flag for enabling large kinematics
  const bool _ld;

private:
  /// method computing the displacement jump and its increment
  void computeJumpInterface();

  /// method computing the interface deformation gradient
  void computeFInterface();

  /// method computing interface rotation and its increment
  void computeRInterface();

  /// method computing deformedinterface normal and its increment
  void computeNormalInterface();

  /// method computing area change and area change increment
  void computeAreaInterface();

  /// method computing the interface strain contribution
  void computeInterfaceStrainRates();

  /// the current and incremental displacement jump
  ///@{
  RealVectorValue _jump;
  RealVectorValue _Djump;
  ///@}

  /// the current and old midplane deformation gradient
  ///@{
  RankTwoTensor _F_average;
  RankTwoTensor _F_average_old;
  ///@}

  /// interface incremental velocity gradient
  RankTwoTensor _DL;

  /// the interface total and incremental rotation
  ///@{
  RankTwoTensor _R_avg;
  RankTwoTensor _DR_avg;
  ///@}

  /// the interface deformed normal and its increment
  ///@{
  RealVectorValue _n_average;
  RealVectorValue _Dn_average;
  ///@}

  /// ratio between deformed and undeformed area
  Real _dadA;
  /// area increment
  Real _Da;
};
