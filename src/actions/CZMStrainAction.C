#include "AddVariableAction.h"
#include "CZMStrainAction.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "RankTwoScalarTools.h"
#include "libmesh/string_to_enum.h"
registerMooseAction("DeerApp", CZMStrainAction, "add_postprocessor");
registerMooseAction("DeerApp", CZMStrainAction, "add_material");
registerMooseAction("DeerApp", CZMStrainAction, "meta_action");

InputParameters
CZMStrainAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Compute the total, normal and sliding strain "
                             "contirbutions of a cohesive model to an RVE.");
  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundary", "The list of boundary names where the cohesive model is present");
  params.addParam<std::vector<SubdomainName>>(
      "block", "The list of subdomain names representing the bulk material");
  params.addParam<PostprocessorName>(
      "bulk_volume_PP", "", "The name of the PP representing the bulk volume at time 0.");
  params.addParam<PostprocessorName>("czm_strain_base_name",
                                     "czm_strain",
                                     "A string containing the base name for interface strain");
  params.addParam<bool>(
      "compute_czm_strain_rate", true, "If true (default) also compute strain rates");
  params.addParam<bool>(
      "compute_equivalent_strain", true, "If true (default) also add all the equivalnet strain");
  MooseEnum strainType("SMALL FINITE", "SMALL");
  params.addParam<MooseEnum>("strain", strainType, "Strain formulation");
  return params;
}

CZMStrainAction::CZMStrainAction(const InputParameters & params)
  : Action(params),
    _block(getParam<std::vector<SubdomainName>>("block")),
    _boundary(getParam<std::vector<BoundaryName>>("boundary")),

    _scaled(!getParam<PostprocessorName>("bulk_volume_PP").empty()),
    _bulk_volume_PP(_scaled ? getParam<PostprocessorName>("bulk_volume_PP") : "czm_strain_V0"),
    _area_ratio_PP("czm_area_ratio"),
    _czm_strain_scale_PP("czm_strain_scale_PP"),
    _czm_strain_base_name(getParam<PostprocessorName>("czm_strain_base_name")),
    _compute_czm_strain_rate(getParam<bool>("compute_czm_strain_rate")),
    _compute_equivalent_strain(getParam<bool>("compute_equivalent_strain"))
{

  // sanity checks
  if (!_scaled && (_block.size() == 0))
    mooseError("CZMStrainAction: The user must provide either bulk_volume_PP "
               "or the list of subdomains representing the bulk material.");

  if (_scaled && (_block.size() > 0))
    mooseError("CZMStrainAction: The user can't provide both the "
               "bulk_volume_PP and the list of subdomains.");
}

void
CZMStrainAction::act()
{
  if (_current_task == "add_postprocessor")
  {
    if (_block.size() > 0)
      computeScalingVolume();
    addInterfaceStrain();
    if (_compute_czm_strain_rate)
      addInterfaceStrainRate();
    if (_compute_equivalent_strain)
    {
      for (unsigned int mp = 0; mp < _czm_mp_strain_names.size(); mp++)
        addEquivalentStrain(_czm_mp_strain_names[mp]);
      if (_compute_czm_strain_rate)
      {
        for (unsigned int mp = 0; mp < _czm_mp_strain_names.size(); mp++)
        {
          auto params_pp = _factory.getValidParams("TimeDerivativePostprocessor");
          params_pp.set<PostprocessorName>("postprocessor") = _czm_mp_strain_names[mp] + "_eq";
          params_pp.set<ExecFlagEnum>("execute_on") = {EXEC_TIMESTEP_END};
          _problem->addPostprocessor(
              "TimeDerivativePostprocessor", _czm_mp_strain_names[mp] + "_eq_rate", params_pp);
        }
      }
    }
  }

  if (_current_task == "add_material")
    addInterfaceStrainMaterial();
}

void
CZMStrainAction::addInterfaceStrain()
{

  for (unsigned int mp = 0; mp < _czm_mp_strain_names.size(); mp++)
    for (auto entry : _tensor_map)
    {
      auto params_pp = _factory.getValidParams("CZMStrainComponent");
      params_pp.set<MaterialPropertyName>("rank_two_tensor") = _czm_mp_strain_names[mp];
      params_pp.set<unsigned int>("index_i") = entry.first.first;
      params_pp.set<unsigned int>("index_j") = entry.first.second;
      params_pp.set<std::vector<BoundaryName>>("boundary") = _boundary;
      params_pp.set<PostprocessorName>("initial_bulk_volume_pp") = _bulk_volume_PP;
      params_pp.set<ExecFlagEnum>("execute_on") = {EXEC_TIMESTEP_END};
      params_pp.set<MooseEnum>("strain") = getParam<MooseEnum>("strain");
      _problem->addPostprocessor(
          "CZMStrainComponent", _czm_mp_strain_names[mp] + "_" + entry.second, params_pp);
    }
}

void
CZMStrainAction::addInterfaceStrainRate()
{
  for (unsigned int mp = 0; mp < _czm_mp_strain_names.size(); mp++)
    for (auto entry : _tensor_map)
    {
      auto params_pp = _factory.getValidParams("TimeDerivativePostprocessor");
      params_pp.set<PostprocessorName>("postprocessor") =
          _czm_mp_strain_names[mp] + "_" + entry.second;
      params_pp.set<ExecFlagEnum>("execute_on") = {EXEC_TIMESTEP_END};
      _problem->addPostprocessor("TimeDerivativePostprocessor",
                                 _czm_mp_strain_names[mp] + "_rate_" + entry.second,
                                 params_pp);
    }
}

void
CZMStrainAction::addEquivalentStrain(const PostprocessorName & rank_two_base_name)
{

  auto params_pp = _factory.getValidParams("RankTwoTensorInvariantPostprocessor");
  params_pp.set<MooseEnum>("invariant") = "EffectiveStrain";
  params_pp.set<PostprocessorName>("rank_two_tensor_base_name") = rank_two_base_name;

  _problem->addPostprocessor(
      "RankTwoTensorInvariantPostprocessor", rank_two_base_name + "_eq", params_pp);
}

void
CZMStrainAction::addInterfaceStrainMaterial()
{
  auto params_mat = _factory.getValidParams("CZMVolumetricStrain");
  params_mat.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params_mat.set<MooseEnum>("strain") = getParam<MooseEnum>("strain");
  _problem->addMaterial("CZMVolumetricStrain", "czm_strains", params_mat);
}

void
CZMStrainAction::computeScalingVolume()
{
  auto params_pp1 = _factory.getValidParams("VolumePostprocessor");
  params_pp1.set<bool>("use_displaced_mesh") = false;
  params_pp1.set<std::vector<SubdomainName>>("block") = _block;
  params_pp1.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;
  _problem->addPostprocessor("VolumePostprocessor", _bulk_volume_PP, params_pp1);
}
