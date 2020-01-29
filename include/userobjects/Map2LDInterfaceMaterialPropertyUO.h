//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Map2LDelem.h"

class Map2LDInterfaceMaterialPropertyUO;

template <> InputParameters validParams<Map2LDInterfaceMaterialPropertyUO>();

/**
 *
 */
class Map2LDInterfaceMaterialPropertyUO : public Map2LDelem {
public:
  static InputParameters validParams();

  Map2LDInterfaceMaterialPropertyUO(const InputParameters &parameters);
  virtual ~Map2LDInterfaceMaterialPropertyUO();

  virtual void initialSetup() override;
  virtual void execute() override;

  Real getValueForLD(dof_id_type ldelem, unsigned int qp) const;
  Real getQpValue(dof_id_type /*elem*/, unsigned int /*side*/,
                  unsigned int /*qp*/) const;

protected:
  virtual Real getMPValue(const unsigned int /*qp*/) const = 0;
  /// this map is used to store QP data.
  std::map<std::pair<dof_id_type, unsigned int>, std::vector<Real>> _map_values;
};
