#include "AddVariableAction.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "RankTwoTensorIntegralAction.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("DeerApp", RankTwoTensorIntegralAction, "add_postprocessor");

InputParameters
RankTwoTensorIntegralAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Compute the volume or area integral for each component of the "
      "provided rank two tensor. Optionally, the compute integrals can be "
      "scaled using an additional postprocessor (scaling_factor_PP), provided "
      "as input.");
  params.addParam<bool>("use_displaced_mesh",
                        true,
                        "If true (default) uses the current volume/area to "
                        "compute the integral.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "rank_two_tensor", "a vector of tensor material property names");
  params.addParam<std::vector<BoundaryName>>(
      "boundary", "The list of boundary names where this action applies");
  params.addParam<std::vector<SubdomainName>>(
      "block",
      "The list of subdomain names where this action applies. If "
      "empty all subdomains are used");
  params.addParam<PostprocessorName>("scaling_factor_PP",
                                     "The name of the PP used as scaling factor.");
  params.addRequiredParam<std::vector<PostprocessorName>>(
      "base_out_names",
      "A vector containing base names of the output "
      "variables, one for each provided rank_two_tensor.");
  return params;
}

RankTwoTensorIntegralAction::RankTwoTensorIntegralAction(const InputParameters & params)
  : Action(params),
    _mp_names(getParam<std::vector<MaterialPropertyName>>("rank_two_tensor")),
    _use_displaced_mesh(getParam<bool>("use_displaced_mesh")),
    _block(params.isParamValid("block") ? getParam<std::vector<SubdomainName>>("block") : std::vector<SubdomainName>()),
    _boundary(params.isParamValid("boundary") ? getParam<std::vector<BoundaryName>>("boundary") : std::vector<BoundaryName>()),
    _scaled(isParamValid("scaling_factor_PP")),
    _scaling_factor_PP(_scaled ? getParam<PostprocessorName>("scaling_factor_PP") : ""),
    _PP_type(!params.isParamValid("boundary")
                 ? (_scaled ? "MaterialTensorIntegralScaled" : "MaterialTensorIntegral")
                 : "MaterialTensorIntegralInterfaceScaled"),
    _base_out_name(getParam<std::vector<PostprocessorName>>("base_out_names"))
{

  /// sanity check
  if (params.isParamValid("block") && params.isParamValid("boundary"))
    mooseError("RankTwoTensorIntegralAction: you can't specify both boundaries "
               "and blocks.");

  if (_base_out_name.size() != _mp_names.size())
    mooseError("RankTwoTensorIntegralAction: The length of rank_two_tensor and "
               "base_out_names input "
               "parameters must be the same! Please check your input file");
}

void
RankTwoTensorIntegralAction::act()
{
  if (_current_task == "add_postprocessor")
  {

    for (unsigned int mp = 0; mp < _mp_names.size(); mp++)
      for (auto entry : _tensor_map)
      {
        auto params_pp = _factory.getValidParams(_PP_type);
        params_pp.set<bool>("use_displaced_mesh") = _use_displaced_mesh;
        params_pp.set<MaterialPropertyName>("rank_two_tensor") = _mp_names[mp];
        params_pp.set<unsigned int>("index_i") = entry.first.first;
        params_pp.set<unsigned int>("index_j") = entry.first.second;

        if (_block.size() > 0)
          params_pp.set<std::vector<SubdomainName>>("block") = _block;
        if (_boundary.size() > 0)
          params_pp.set<std::vector<BoundaryName>>("boundary") = _boundary;
        if (_scaled)
          params_pp.set<PostprocessorName>("scaling_factor_PP") = _scaling_factor_PP;

        _problem->addPostprocessor(_PP_type, _base_out_name[mp] + "_" + entry.second, params_pp);
      }
  }
}
