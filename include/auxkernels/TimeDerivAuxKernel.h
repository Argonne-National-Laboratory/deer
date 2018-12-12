#ifndef TIMEDERIVAUXKERNEL_H
#define TIMEDERIVAUXKERNEL_H

#include "AuxKernel.h"

class TimeDerivAuxKernel;

template <>
InputParameters validParams<TimeDerivAuxKernel>();

class TimeDerivAuxKernel: public AuxKernel
{
 public:
  TimeDerivAuxKernel(const InputParameters & parameters);

 protected:
  virtual Real computeValue() override;

 protected:
  const VariableValue & _coupled_new;
  const VariableValue & _coupled_old;
};


#endif // TIMEDERIVAUXKERNEL_H
