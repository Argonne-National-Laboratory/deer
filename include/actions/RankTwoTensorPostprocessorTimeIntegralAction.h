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

  // map between tensor components and names
  const std::map<std::pair<int, int>, std::string> _tensor_map = {
      {std::make_pair(0, 0), "x-x"}, {std::make_pair(1, 1), "y-y"},
      {std::make_pair(2, 2), "z-z"}, {std::make_pair(0, 1), "x-y"},
      {std::make_pair(0, 2), "x-z"}, {std::make_pair(1, 2), "y-z"}};
};
