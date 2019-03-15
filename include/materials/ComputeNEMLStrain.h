#ifndef COMPUTENEMLSTRAIN_H
#define COMPUTENEMLSTRAIN_H

#include "RankTwoTensor.h"
#include "ComputeNEMLStrainBase.h"

class ComputeNEMLStrain;

template <>
InputParameters validParams<ComputeNEMLStrain>();

class ComputeNEMLStrain: public ComputeNEMLStrainBase
{
 public:
  ComputeNEMLStrain(const InputParameters & parameters);
  virtual ~ComputeNEMLStrain() {};

 protected:
  virtual void computeQpProperties() override;
};

#endif // COMPUTENEMLSTRAIN_H
