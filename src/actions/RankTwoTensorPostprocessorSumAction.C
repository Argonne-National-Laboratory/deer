#include "FEProblem.h"
#include "MooseMesh.h"
#include "RankTwoTensorPostprocessorSumAction.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("DeerApp", RankTwoTensorPostprocessorSumAction,
                    "add_postprocessor");

InputParameters RankTwoTensorPostprocessorSumAction::validParams() {
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Compute the tensor sum of pairs of postprocessors "
      "each representing a component of a rank two tensor");
  params.addRequiredParam<std::vector<PostprocessorName>>(
      "pp_base_names_1", "a vector of postprocessor base names represting the "
                         "first value of the sum.");
  params.addRequiredParam<std::vector<PostprocessorName>>(
      "pp_base_names_2", "a vector of postprocessor base names represting the "
                         "second value of the sum.");
  params.addRequiredParam<std::vector<PostprocessorName>>(
      "base_out_names",
      "A vector containing base names of the output "
      "variables, one for each provided postprocessor base name.");
  return params;
}

RankTwoTensorPostprocessorSumAction::RankTwoTensorPostprocessorSumAction(
    const InputParameters &params)
    : Action(params), _pp_base_names_1(getParam<std::vector<PostprocessorName>>(
                          "pp_base_names_1")),
      _pp_base_names_2(
          getParam<std::vector<PostprocessorName>>("pp_base_names_2")),
      _base_out_name(
          getParam<std::vector<PostprocessorName>>("base_out_names")) {

  /// sanity check
  if (_base_out_name.size() != _pp_base_names_1.size())
    mooseError("RankTwoTensorPostprocessorSumAction: The length of "
               "pp_base_names_1 and "
               "base_out_names input "
               "parameters must be the same! Please check your input file");

  if (_base_out_name.size() != _pp_base_names_2.size())
    mooseError("RankTwoTensorPostprocessorSumAction: The length of "
               "pp_base_names_2 and "
               "base_out_names input "
               "parameters must be the same! Please check your input file");

  if (_pp_base_names_1.size() != _pp_base_names_2.size())
    mooseError("RankTwoTensorPostprocessorSumAction: The length of "
               "pp_base_names_1 and "
               "pp_base_names_2 input "
               "parameters must be the same! Please check your input file");
}

void RankTwoTensorPostprocessorSumAction::act() {
  if (_current_task == "add_postprocessor") {

    for (unsigned int pp = 0; pp < _pp_base_names_1.size(); pp++)
      for (auto entry : tensor_map) {
        auto params_pp = _factory.getValidParams("SumPostprocessor");
        params_pp.set<PostprocessorName>("value1") =
            _pp_base_names_1[pp] + "_" + entry.second;
        params_pp.set<PostprocessorName>("value2") =
            _pp_base_names_2[pp] + "_" + entry.second;

        _problem->addPostprocessor("SumPostprocessor",
                                   _base_out_name[pp] + "_" + entry.second,
                                   params_pp);
      }
  }
}
