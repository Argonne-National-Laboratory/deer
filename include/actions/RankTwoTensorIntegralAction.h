#pragma once

#include "Action.h"

class RankTwoTensorIntegralAction : public Action {
public:
  static InputParameters validParams();
  RankTwoTensorIntegralAction(const InputParameters &params);

  virtual void act();

protected:
  const std::vector<MaterialPropertyName> _mp_names;
  const bool _use_displaced_mesh;
  const std::vector<SubdomainName> _block;
  const std::vector<BoundaryName> _boundary;
  const bool _scaled;
  const PostprocessorName _scaling_factor_PP;
  const std::string _PP_type;
  const std::vector<PostprocessorName> _base_out_name;

  // map between tensor components and names
  const std::map<std::pair<int, int>, std::string> _tensor_map = {
      {std::make_pair(0, 0), "xx"}, {std::make_pair(1, 1), "yy"},
      {std::make_pair(2, 2), "zz"}, {std::make_pair(0, 1), "xy"},
      {std::make_pair(0, 2), "xz"}, {std::make_pair(1, 2), "yz"}};
};
