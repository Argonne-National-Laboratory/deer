//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"

class LatticeTransformGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  LatticeTransformGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

private:
  void setupTransforms();

protected:
  /// the three lattice vectors
  Point _a1, _a2, _a3;

  /// the input mesh
  std::unique_ptr<MeshBase> & _input;

  /// Transformation matrix
  RankTwoTensor _T;

  /// Shift
  Point _b;
};
