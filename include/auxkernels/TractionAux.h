//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MaterialAuxBase.h"
#include "MooseEnum.h"
#include "RankTwoTensor.h"
class TractionAux;

/**
 * TractionAux is designed to take the data in the RankTwoTensor
 * material property, for example stress or strain, and output the value for the
 * supplied indices.
 */
class TractionAux : public MaterialAuxBase<RankTwoTensor> {
public:
  static InputParameters validParams();
  TractionAux(const InputParameters &parameters);

protected:
  virtual void initialSetup() override;
  virtual Real getRealValue() override;
  const MooseArray<Point> &_normals;
  const bool _PK1;
  const bool _large_kinematics;
  /// number of displacement components
  const unsigned int _ndisp;

  /// the coupled displacement gradient
  std::vector<const VariableGradient *> _grad_disp;

  enum class scalarType {
    normal,
    shear1,
    shear2,
    shear_norm,
    X,
    Y,
    Z
  } _scalar_type;
};
