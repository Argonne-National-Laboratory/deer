#pragma once

#include "Action.h"

class RankTwoTensorPostprocessorSumAction : public Action {
public:
  static InputParameters validParams();
  RankTwoTensorPostprocessorSumAction(const InputParameters &params);

  virtual void act();

protected:
  const std::vector<PostprocessorName> _pp_base_names_1;
  const std::vector<PostprocessorName> _pp_base_names_2;
  const std::vector<PostprocessorName> _base_out_name;
};
