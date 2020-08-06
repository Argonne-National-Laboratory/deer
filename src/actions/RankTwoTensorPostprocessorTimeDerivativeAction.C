#include "FEProblem.h"
#include "MooseMesh.h"
#include "RankTwoTensorPostprocessorTimeDerivativeAction.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("DeerApp", RankTwoTensorPostprocessorTimeDerivativeAction,
                    "add_postprocessor");

InputParameters RankTwoTensorPostprocessorTimeDerivativeAction::validParams() {
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

RankTwoTensorPostprocessorTimeDerivativeAction::
    RankTwoTensorPostprocessorTimeDerivativeAction(
        const InputParameters &params)
    : Action(params),
      _pp_base_names(getParam<std::vector<PostprocessorName>>("pp_base_names")),
      _base_out_name(
          getParam<std::vector<PostprocessorName>>("base_out_names")) {

  /// sanity check
  if (_base_out_name.size() != _pp_base_names.size())
    mooseError("RankTwoTensorPostprocessorTimeDerivativeAction: The length of "
               "pp_base_names and "
               "base_out_names input "
               "parameters must be the same! Please check your input file");
}

void RankTwoTensorPostprocessorTimeDerivativeAction::act() {
  if (_current_task == "add_postprocessor") {

    for (unsigned int pp = 0; pp < _pp_base_names.size(); pp++)
      for (auto entry : tensor_map) {
        auto params_pp = _factory.getValidParams("TimeDerivativePostprocessor");
        params_pp.set<PostprocessorName>("value") =
            _pp_base_names[pp] + "_" + entry.second;

        _problem->addPostprocessor("TimeDerivativePostprocessor",
                                   _base_out_name[pp] + "_" + entry.second,
                                   params_pp);
      }
  }
}
