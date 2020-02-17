#pragma once

#include "ComputeNEMLStrainBase.h"
#include "RankTwoTensor.h"

class ComputeNEMLStrain : public ComputeNEMLStrainBase {
public:
  static InputParameters validParams();
  ComputeNEMLStrain(const InputParameters &parameters);
  virtual ~ComputeNEMLStrain(){};

protected:
  virtual void computeQpProperties() override;
};
