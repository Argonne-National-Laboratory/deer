#pragma once

#include "Kernel.h"

class AdiabaticHeating : public Kernel
{
public:
  static InputParameters validParams();
  AdiabaticHeating(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

protected:
  const VariableValue & _heat;
  const Real _fraction;
};
