//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

class ScalarPrincipalValues : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  ScalarPrincipalValues(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual PostprocessorValue getValue() const;

protected:
  const VariableValue & _scalar_var;
  const Order & _scalar_order;
  const size_t _rank;
  PostprocessorValue _value;
};
