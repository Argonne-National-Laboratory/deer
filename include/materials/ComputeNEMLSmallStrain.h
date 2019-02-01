#ifndef COMPUTENEMLSMALLSTRAIN_H
#define COMPUTENEMLSMALLSTRAIN_H

#include "RankTwoTensor.h"
#include "ComputeNEMLStrainBase.h"

class ComputeNEMLSmallStrain;

template <>
InputParameters validParams<ComputeNEMLSmallStrain>();

class ComputeNEMLSmallStrain: public ComputeNEMLStrainBase
{
 public:
  ComputeNEMLSmallStrain(const InputParameters & parameters);
  virtual ~ComputeNEMLSmallStrain() {};

 protected:
  virtual void computeProperties() override;
};

#endif // COMPUTENEMLSMALLSTRAIN_H
