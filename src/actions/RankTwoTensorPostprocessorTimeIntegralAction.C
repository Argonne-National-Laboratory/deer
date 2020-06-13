#include "FEProblem.h"
#include "MooseMesh.h"
#include "RankTwoTensorPostprocessorTimeIntegralAction.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("DeerApp", RankTwoTensorPostprocessorTimeIntegralAction,
                    "add_postprocessor");

InputParameters RankTwoTensorPostprocessorTimeIntegralAction::validParams() {
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Compute the time integral for each component of several postprocessors "
      "each representing a component of a rank two tensor");
  params.addRequiredParam<std::vector<PostprocessorName>>(
      "pp_base_names", "a vector of postprocessor base names.");
  params.addRequiredParam<std::vector<PostprocessorName>>(
      "base_out_names",
      "A vector containing base names of the output "
      "variables, one for each provided postprocessor base name.");
  return params;
}

RankTwoTensorPostprocessorTimeIntegralAction::
    RankTwoTensorPostprocessorTimeIntegralAction(const InputParameters &params)
    : Action(params),
      _pp_base_names(getParam<std::vector<PostprocessorName>>("pp_base_names")),
      _base_out_name(
          getParam<std::vector<PostprocessorName>>("base_out_names")) {

  /// sanity check
  if (_base_out_name.size() != _pp_base_names.size())
    mooseError("RankTwoTensorPostprocessorTimeIntegralAction: The length of "
               "pp_base_names and "
               "base_out_names input "
               "paramters must be the same! Please check your input file");
}

void RankTwoTensorPostprocessorTimeIntegralAction::act() {
  if (_current_task == "add_postprocessor") {

    for (unsigned int pp = 0; pp < _pp_base_names.size(); pp++)
      for (auto entry : tensor_map) {
        auto params_pp = _factory.getValidParams("TimeIntegralPostprocessor");
        params_pp.set<PostprocessorName>("postprocessor") =
            _pp_base_names[pp] + "_" + entry.second;

        _problem->addPostprocessor("TimeIntegralPostprocessor",
                                   _base_out_name[pp] + "_" + entry.second,
                                   params_pp);
      }
  }
}
