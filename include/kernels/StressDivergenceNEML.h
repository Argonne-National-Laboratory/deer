#ifndef STRESSDIVERGENCENEML_H
#define STRESSDIVERGENCENEML_H

#include "Kernel.h"

class StressDivergenceNEML;

template <>
InputParameters validParams<StressDivergenceNEML>();

class StressDivergenceNEML: public DerivativeMaterialInterface<Kernel>
{
 public:
  StressDivergenceNEML(const InputParameters & parameters);
  virtual ~StressDivergenceNEML() {};
  
 protected:
  // These must be overwritten for Bbar-type stabilizations
  virtual void precalculateResidual() override;
  virtual void precalculateJacobian() override;
  virtual void precalculateOffDiagJacobian(unsigned int jvar);

  // These are the standard implementation functions
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
};

#endif // STRESSDIVERGENCENEML_H
