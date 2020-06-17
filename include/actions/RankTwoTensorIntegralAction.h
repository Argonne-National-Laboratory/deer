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
};
