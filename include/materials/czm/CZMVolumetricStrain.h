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
 * This is the base Material class for implementing a traction separation
 * material model including large rotation and area changes. The responsibility
 * of this class is to rotate the displacement jump from global to local
 * coordinate and rotate back traction and traction derivatives. The local
 * coordinate system assumes the following order: opening, tangential1,
 * tangential2. Note that tangential1, tangential2 are arbitrary and therefore
 * the interface assumes an in-plane isotropic behavior. By overriding
 * computeTractionIncrementAndDerivatives in aderived class, different traction
 * separation laws can be implemented. The
 * computeTractionIncrementAndDerivatives method assumes calculations are
 * performed in the local frame. CZM laws should always be implemented in 3D
 * even if they are going to be used in 2D or 1D simulations. This class assumes
 * large rotation, large area changes, and that the traction separation law is
 * only dependent upon the the displacement jump.
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

  /// the value of the traction in global and local coordinates
  ///@{
  MaterialProperty<RankTwoTensor> &_czm_total_strain_rate;
  MaterialProperty<RankTwoTensor> &_czm_normal_strain_rate;
  MaterialProperty<RankTwoTensor> &_czm_sliding_strain_rate;
  ///@}

  const bool _large_kinematics;
};
