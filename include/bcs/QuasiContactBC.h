#pragma once

#include "IntegratedBC.h"
#include "Function.h"

class QuasiContactBC: public IntegratedBC
{
 public:
  static InputParameters validParams();

  QuasiContactBC(const InputParameters & parameters);

 protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

 private:
  const Function & _func;
  const Real & _v;
};
