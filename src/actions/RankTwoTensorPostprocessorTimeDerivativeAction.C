#include "FEProblem.h"
#include "MooseMesh.h"
#include "RankTwoTensorPostprocessorTimeDerivativeAction.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("DeerApp", RankTwoTensorPostprocessorTimeDerivativeAction, "add_postprocessor");

InputParameters
RankTwoTensorPostprocessorTimeDerivativeAction::validParams()
{
  InputParameters params = RankTwoTensorPostprocessorTimeIntegralAction::validParams();
  params.addClassDescription("Compute the time derivative for each component of several "
                             "postprocessors  each representing a component of a rank two tensor");
  return params;
}

RankTwoTensorPostprocessorTimeDerivativeAction::RankTwoTensorPostprocessorTimeDerivativeAction(
    const InputParameters & params)
  : RankTwoTensorPostprocessorTimeIntegralAction(params)
{
}

void
RankTwoTensorPostprocessorTimeDerivativeAction::act()
{
  if (_current_task == "add_postprocessor")
  {

    for (unsigned int pp = 0; pp < _pp_base_names.size(); pp++)
      for (auto entry : _tensor_map)
      {
        auto params_pp = _factory.getValidParams("TimeDerivativePostprocessor");
        params_pp.set<PostprocessorName>("postprocessor") = _pp_base_names[pp] + "_" + entry.second;

        _problem->addPostprocessor(
            "TimeDerivativePostprocessor", _base_out_name[pp] + "_" + entry.second, params_pp);
      }
  }
}
