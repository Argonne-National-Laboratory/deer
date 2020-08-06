#include "AddVariableAction.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "RankTwoTensorIntegralAction.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("DeerApp", RankTwoTensorIntegralAction,
                    "add_postprocessor");

InputParameters RankTwoTensorIntegralAction::validParams() {
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Compute the volume or area integral for each component of the "
      "provided rank two tensor. Optionally, the compute integrals can be "
      "scaled using an additional postprocessor (scaling_factor_PP), provided "
      "as input.");
  params.addParam<bool>("use_displaced_mesh", true,
                        "If true (default) uses the current volume/area to "
                        "compute the integral.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "rank_two_tensor", "a vector of tensor material property names");
  params.addParam<std::vector<BoundaryName>>(
      "boundary", "The list of boundary names where this action applies");
  params.addParam<std::vector<SubdomainName>>(
      "block", "The list of subdomain names where this action applies. If "
               "empty all subdomains are used");
  params.addParam<PostprocessorName>(
      "scaling_factor_PP", "The name of the PP used as scaling factor.");
  params.addRequiredParam<std::vector<PostprocessorName>>(
      "base_out_names", "A vector containing base names of the output "
                        "variables, one for each provided rank_two_tensor.");
  params.addParam<bool>("normalize_integral_by_area", false,
                        "If true normalize the integral by the interface area. "
                        "This is done in addition to the scaling_factor_PP.");
  params.addParam<bool>("czm", false,
                        "If true use integration and area normalization are "
                        "performed over the true cohesive zone area. Should be "
                        "true when using cohesive zone wiuth alrge kinematics");
  return params;
}

RankTwoTensorIntegralAction::RankTwoTensorIntegralAction(
    const InputParameters &params)
    : Action(params),
      _mp_names(getParam<std::vector<MaterialPropertyName>>("rank_two_tensor")),
      _use_displaced_mesh(getParam<bool>("use_displaced_mesh")),
      _block(getParam<std::vector<SubdomainName>>("block")),
      _boundary(getParam<std::vector<BoundaryName>>("boundary")),

      _scaled(isParamValid("scaling_factor_PP")),
      _scaling_factor_PP(
          _scaled ? getParam<PostprocessorName>("scaling_factor_PP") : ""),
      _czm(getParam<bool>("czm")),
      _PP_type(_boundary.size() == 0
                   ? (_scaled ? "MaterialTensorIntegralScaled"
                              : "MaterialTensorIntegral")
                   : _czm ? "MaterialTensorIntegralCZMScaled"
                          : "MaterialTensorIntegralInterfaceScaled"),
      _base_out_name(
          getParam<std::vector<PostprocessorName>>("base_out_names")),
      _normalize_integral_by_area(
          getParam<bool>("normalize_integral_by_area")) {

  /// sanity check
  if (params.isParamValid("block") && params.isParamValid("boundary"))
    mooseError("RankTwoTensorIntegralAction: you can't specify both boundaries "
               "and blocks.");

  if (_base_out_name.size() != _mp_names.size())
    mooseError("RankTwoTensorIntegralAction: The length of rank_two_tensor and "
               "base_out_names input "
               "parameters must be the same! Please check your input file");

  if (!params.isParamValid("boundary") && _czm)
    mooseError(
        "If czm paramter is true than the paramter boundary cannot be empty");
}

void RankTwoTensorIntegralAction::act() {
  if (_current_task == "add_postprocessor") {

    for (unsigned int mp = 0; mp < _mp_names.size(); mp++)
      for (auto entry : tensor_map) {
        auto params_pp = _factory.getValidParams(_PP_type);
        params_pp.set<bool>("use_displaced_mesh") = _use_displaced_mesh;
        params_pp.set<MaterialPropertyName>("rank_two_tensor") = _mp_names[mp];
        params_pp.set<unsigned int>("index_i") = entry.first.first;
        params_pp.set<unsigned int>("index_j") = entry.first.second;

        if (_block.size() > 0)
          params_pp.set<std::vector<SubdomainName>>("block") = _block;
        if (_boundary.size() > 0) {
          params_pp.set<std::vector<BoundaryName>>("boundary") = _boundary;
          params_pp.set<bool>("normalize_integral_by_area") =
              _normalize_integral_by_area;
        }
        if (_scaled)
          params_pp.set<PostprocessorName>("scaling_factor_PP") =
              _scaling_factor_PP;

        _problem->addPostprocessor(
            _PP_type, _base_out_name[mp] + "_" + entry.second, params_pp);
      }
  }
}
