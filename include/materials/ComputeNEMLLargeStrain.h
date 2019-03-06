#ifndef COMPUTENEMLLARGESTRAIN_H
#define COMPUTENEMLLARGESTRAIN_H

#include "RankTwoTensor.h"
#include "ComputeNEMLStrainBase.h"

class ComputeNEMLLargeStrain;

template <>
InputParameters validParams<ComputeNEMLLargeStrain>();

class ComputeNEMLLargeStrain: public ComputeNEMLStrainBase
{
 public:
  ComputeNEMLLargeStrain(const InputParameters & parameters);
  virtual ~ComputeNEMLLargeStrain() {};

 protected:
  virtual void computeQpProperties() override;
};

#endif // COMPUTENEMLLARGESTRAIN_H
