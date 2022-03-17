//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 * Empirical reaction heating 
 */
class EmpiricalReactionHeat : public Material {
 public:
  static InputParameters validParams();
  EmpiricalReactionHeat(const InputParameters &parameters);
 
 protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;
 
 protected:
  Real _W0, _T0, _tau;

  MaterialProperty<bool> &_active;
  const MaterialProperty<bool> &_active_old;

  MaterialProperty<Real> &_power;
  const MaterialProperty<Real> & _power_old;

  const VariableValue & _temperature;
};
