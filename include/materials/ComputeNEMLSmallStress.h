#ifndef COMPUTENEMLSMALLSTRESS_H
#define COMPUTENEMLSMALLSTRESS_H

#include "ComputeNEMLStressBase.h"

class ComputeNEMLSmallStress;

template <>
InputParameters validParams<ComputeNEMLSmallStress>();

class ComputeNEMLSmallStress: public ComputeNEMLStressBase
{
 public:
  ComputeNEMLSmallStress(const InputParameters & parameters);
  virtual ~ComputeNEMLSmallStress() {};

 protected:
  virtual void stressUpdate();

 private:
  void updateStrain();
};

#endif
