//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceKernel.h"

/// DG kernel implementing cohesive zone models (CZM) for a 1D/2D/3D traction
/// separation laws based on the displacement jump. This kernel operates only on
/// a single displacement compenent.
/// One kernel is required for each mesh dimension.
class CZMInterfaceKernelPK1 : public InterfaceKernel {
public:
  static InputParameters validParams();
  CZMInterfaceKernelPK1(const InputParameters &parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type,
                                        unsigned int jvar);

  /// method computing the jacobian contribution due to rotations and area
  /// changes
  Real JacLD(const unsigned int cc, const bool neighbor) const;

  /// the displacement component this kernel is operating on (0=x, 1=y, 2 =z)
  const unsigned int _component;

  /// number of displacement components
  const unsigned int _ndisp;

  /// Coupled displacement component variable IDs
  ///@{
  std::vector<unsigned int> _disp_var;
  std::vector<unsigned int> _disp_neighbor_var;
  ///@}

  // pointer to displacement variables
  std::vector<MooseVariable *> _vars;

  // values of the traction and traction derivatives used by the kernel
  ///@{
  const MaterialProperty<RealVectorValue> &_PK1traction;
  const MaterialProperty<RankTwoTensor> &_dPK1traction_djumpglobal;
  const MaterialProperty<RankThreeTensor> &_dPK1traction_dF;
  ///@}
};
