#pragma once

#include "RankTwoTensorPostprocessorTimeIntegralAction.h"

class RankTwoTensorPostprocessorTimeDerivativeAction
  : public RankTwoTensorPostprocessorTimeIntegralAction
{
public:
  static InputParameters validParams();
  RankTwoTensorPostprocessorTimeDerivativeAction(const InputParameters & params);

  virtual void act() override;
};
