#ifndef TIMEDERIVAUXKERNEL_H
#define TIMEDERIVAUXKERNEL_H

#include "AuxKernel.h"

class TimeDerivAuxKernel;

class TimeDerivAuxKernel : public AuxKernel {
public:
  static InputParameters validParams();
  TimeDerivAuxKernel(const InputParameters &parameters);

protected:
  virtual Real computeValue() override;

protected:
  const VariableValue &_coupled_new;
  const VariableValue &_coupled_old;
};

#endif // TIMEDERIVAUXKERNEL_H
