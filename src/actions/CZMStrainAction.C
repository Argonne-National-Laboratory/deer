#include "AddVariableAction.h"
#include "CZMStrainAction.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("DeerApp", CZMStrainAction, "add_postprocessor");
registerMooseAction("DeerApp", CZMStrainAction, "add_material");
registerMooseAction("DeerApp", CZMStrainAction, "meta_action");

InputParameters CZMStrainAction::validParams() {
  InputParameters params = Action::validParams();
  params.addClassDescription("Compute the total, normal and sliding strain "
                             "contirbutions of a cohesive model to an RVE.");
  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundary",
      "The list of boundary names where the cohesive model is present");
  params.addParam<std::vector<SubdomainName>>(
      "block", "The list of subdomain names representing the bulk material");
  params.addParam<PostprocessorName>(
      "bulk_volume_PP", "",
      "The name of the PP representing the bulk volume at time 0.");
  params.addParam<PostprocessorName>(
      "czm_strain_base_name", "czm_strain",
      "A string containing the base name for interface strain");
  params.addParam<bool>(
      "add_integrated_interface_strains", true,
      "If true (default) also add all integrated interface strains");
  return params;
}

CZMStrainAction::CZMStrainAction(const InputParameters &params)
    : Action(params),

      _block(getParam<std::vector<SubdomainName>>("block")),
      _boundary(getParam<std::vector<BoundaryName>>("boundary")),

      _scaled(!getParam<PostprocessorName>("bulk_volume_PP").empty()),
      _bulk_volume_PP(_scaled ? getParam<PostprocessorName>("bulk_volume_PP")
                              : "czm_strain_V0"),

      _czm_strain_base_name(
          getParam<PostprocessorName>("czm_strain_base_name")),
      _add_integrated_interface_strains(
          getParam<bool>("add_integrated_interface_strains")) {

  /// sanity checks
  if (!_scaled && (_block.size() == 0))
    mooseError("CZMStrainAction: The user must provide either bulk_volume_PP "
               "or the list of subdomains representing the bulk material.");

  if (_scaled && (_block.size() > 0))
    mooseError("CZMStrainAction: The user can't provide both the "
               "bulk_volume_PP and the list of subdomains.");
}

void CZMStrainAction::act() {
  if (_current_task == "add_postprocessor") {
    if (_block.size() > 0) {
      auto params_pp = _factory.getValidParams("VolumePostprocessor");
      params_pp.set<bool>("use_displaced_mesh") = false;
      params_pp.set<std::vector<SubdomainName>>("block") = _block;
      params_pp.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;
      _problem->addPostprocessor("VolumePostprocessor", _bulk_volume_PP,
                                 params_pp);
    }
  }

  if (_current_task == "meta_action") {
    addInterfaceStrainRateAction();
    if (_add_integrated_interface_strains)
      addIntegrateInterfaceStrainRateAction();
  }
}

void CZMStrainAction::addInterfaceStrainRateAction() {
  auto params_ma =
      _action_factory.getValidParams("RankTwoTensorIntegralAction");
  params_ma.set<ActionWarehouse *>("awh") = &_awh;

  params_ma.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params_ma.set<bool>("use_displaced_mesh") = false;
  params_ma.set<std::vector<MaterialPropertyName>>("rank_two_tensor") = {
      "czm_total_strain_rate", "czm_normal_strain_rate",
      "czm_sliding_strain_rate"};
  params_ma.set<std::vector<PostprocessorName>>("base_out_names") = {
      _czm_strain_base_name + "_total_rate",
      _czm_strain_base_name + "_normal_rate",
      _czm_strain_base_name + "_sliding_rate"};
  params_ma.set<PostprocessorName>("scaling_factor_PP") = _bulk_volume_PP;

  auto action = _action_factory.create("RankTwoTensorIntegralAction",
                                       "RankTwoTensorIntegralAction/czm_strain",
                                       params_ma);

  _awh.addActionBlock(action);
}

void CZMStrainAction::addIntegrateInterfaceStrainRateAction() {
  auto params_ti = _action_factory.getValidParams(
      "RankTwoTensorPostprocessorTimeIntegralAction");
  params_ti.set<ActionWarehouse *>("awh") = &_awh;
  params_ti.set<std::vector<PostprocessorName>>("pp_base_names") = {
      _czm_strain_base_name + "_total_rate",
      _czm_strain_base_name + "_normal_rate",
      _czm_strain_base_name + "_sliding_rate"};
  params_ti.set<std::vector<PostprocessorName>>("base_out_names") = {
      _czm_strain_base_name + "_total", _czm_strain_base_name + "_normal",
      _czm_strain_base_name + "_sliding"};

  auto action_ti = _action_factory.create(
      "RankTwoTensorPostprocessorTimeIntegralAction",
      "RankTwoTensorPostprocessorTimeIntegralAction/czm_strain", params_ti);

  _awh.addActionBlock(action_ti);
}
