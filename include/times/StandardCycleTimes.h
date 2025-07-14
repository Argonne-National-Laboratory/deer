//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Moose includes
#include "Times.h"

/**
 * Simple times from an input parameter
 */
class StandardCycleTimes : public Times
{
public:
  static InputParameters validParams();
  StandardCycleTimes(const InputParameters & parameters);
  virtual ~StandardCycleTimes() = default;

protected:
  virtual void initialize() override {}

protected:
  Real _load, _tension_hold, _compression_hold;
  size_t _number_load, _number_tension_hold, _number_compression_hold, _cycles;
};
