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

template <> InputParameters validParams<TractionAux>();

class TractionAux : public MaterialAuxBase<RankTwoTensor> {
public:
  TractionAux(const InputParameters &parameters);

protected:
  virtual Real getRealValue() override;
  const MooseArray<Point> &_normals;
  const MooseEnum _scalar_type;

private:
  MooseEnum traction_component();
};
