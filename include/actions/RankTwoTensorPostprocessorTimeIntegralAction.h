#pragma once

#include "Action.h"

class RankTwoTensorPostprocessorTimeIntegralAction : public Action {
public:
  static InputParameters validParams();
  RankTwoTensorPostprocessorTimeIntegralAction(const InputParameters &params);

  virtual void act();

protected:
  const std::vector<PostprocessorName> _pp_base_names;
  const std::vector<PostprocessorName> _base_out_name;
};
