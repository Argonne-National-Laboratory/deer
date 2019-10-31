#include "AddVariableAction.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "NEMLMechanicsAction.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("DeerApp", NEMLMechanicsAction, "add_variable");
registerMooseAction("DeerApp", NEMLMechanicsAction, "add_kernel");
registerMooseAction("DeerApp", NEMLMechanicsAction, "add_material");
registerMooseAction("DeerApp", NEMLMechanicsAction, "add_aux_variable");
registerMooseAction("DeerApp", NEMLMechanicsAction, "add_aux_kernel");

const std::vector<std::string> all_tensors = {
    "mechanical_strain", "stress", "elastic_strain", "inelastic_strain"};
const std::vector<std::string> all_scalars = {"energy", "dissipation"};

const std::map<std::pair<int, int>, std::string> tensor_map = {
    {std::make_pair(0, 0), "x-x"}, {std::make_pair(1, 1), "y-y"},
    {std::make_pair(2, 2), "z-z"}, {std::make_pair(0, 1), "x-y"},
    {std::make_pair(0, 2), "x-z"}, {std::make_pair(1, 2), "y-z"}};

template <> InputParameters validParams<NEMLMechanicsAction>() {
  InputParameters params = validParams<Action>();

  params.addRequiredParam<std::vector<VariableName>>(
      "displacements", "The displacement variables");
  params.addParam<bool>("add_displacements", true,
                        "Add the displacement variables");

  MooseEnum kinematicType("small large", "small");
  params.addParam<MooseEnum>("kinematics", kinematicType,
                             "Kinematic formulation");

  params.addParam<bool>(
      "add_all_output", false,
      "Dump all the usual stress and strain variables to the output");
  params.addParam<bool>(
      "add_all_output", false,
      "Dump all the usual stress and strain variables to the output");

  params.addParam<std::vector<MaterialPropertyName>>("eigenstrains",
                                                     std::vector<MaterialPropertyName>(),
                                                     "Names of the eigenstrains");

  return params;
}

NEMLMechanicsAction::NEMLMechanicsAction(const InputParameters &params)
    : Action(params),
      _displacements(getParam<std::vector<VariableName>>("displacements")),
      _ndisp(_displacements.size()),
      _add_disp(getParam<bool>("add_displacements")),
      _add_all(getParam<bool>("add_all_output")),
      _kinematics(getParam<MooseEnum>("kinematics").getEnum<Kinematics>()),
      _eigenstrains(getParam<std::vector<MaterialPropertyName>>("eigenstrains"))
{

}

void NEMLMechanicsAction::act() {
  if (_current_task == "add_variable") {
    const bool second = _problem->mesh().hasSecondOrderElements();

    InputParameters params = _factory.getValidParams("MooseVariableBase");
    params.set<MooseEnum>("family") = "LAGRANGE";
    if (second)
      params.set<MooseEnum>("order") = "SECOND";
    else
      params.set<MooseEnum>("order") = "FIRST";

    for (const auto &disp_name : _displacements) {
      auto fe_type = AddVariableAction::feType(params);
      auto var_type = AddVariableAction::determineType(fe_type, 1);
      _problem->addVariable(var_type, disp_name, params);
    }
  } else if (_current_task == "add_material") {
    // Add the strain calculator
    auto params = _factory.getValidParams("ComputeNEMLStrain");

    params.set<std::vector<VariableName>>("displacements") = _displacements;
    params.set<std::vector<MaterialPropertyName>>("eigenstrain_names") = _eigenstrains;
    params.set<bool>("large_kinematics") = _kin_mapper[_kinematics];

    _problem->addMaterial("ComputeNEMLStrain", "strain", params);
  } else if (_current_task == "add_kernel") {
    // Add the kernels
    for (unsigned int i = 0; i < _ndisp; ++i) {
      auto params = _factory.getValidParams("StressDivergenceNEML");

      params.set<std::vector<VariableName>>("displacements") = _displacements;
      params.set<NonlinearVariableName>("variable") = _displacements[i];
      params.set<unsigned int>("component") = i;
      params.set<bool>("use_displaced_mesh") = _kin_mapper[_kinematics];

      std::string name = "SD_" + Moose::stringify(i);

      _problem->addKernel("StressDivergenceNEML", name, params);
    }
  } else if (_current_task == "add_aux_variable") {
    if (_add_all) {
      for (auto name : all_tensors)
        _add_tensor_variable(name);
      for (auto name : all_scalars)
        _add_scalar_variable(name);
    }
  } else if (_current_task == "add_aux_kernel") {
    if (_add_all) {
      for (auto name : all_tensors)
        _add_tensor_aux(name);
      for (auto name : all_scalars)
        _add_scalar_aux(name);
    }
  }
}

void NEMLMechanicsAction::_add_tensor_variable(std::string name) {
  for (auto entry : tensor_map) {
    _add_scalar_variable(name + "_" + entry.second);
  }
}

void NEMLMechanicsAction::_add_scalar_variable(std::string name) {

  InputParameters params = _factory.getValidParams("MooseVariableBase");
  params.set<MooseEnum>("family") = "MONOMIAL";
  params.set<MooseEnum>("order") = "CONSTANT";

  auto fe_type = AddVariableAction::feType(params);
  auto var_type = AddVariableAction::determineType(fe_type, 1);

  _problem->addAuxVariable(var_type, name, params);
}

void NEMLMechanicsAction::_add_tensor_aux(std::string name) {
  for (auto entry : tensor_map) {
    auto params = _factory.getValidParams("RankTwoAux");

    params.set<MaterialPropertyName>("rank_two_tensor") = name;
    params.set<AuxVariableName>("variable") = name + "_" + entry.second;
    params.set<unsigned int>("index_i") = entry.first.first;
    params.set<unsigned int>("index_j") = entry.first.second;

    _problem->addAuxKernel("RankTwoAux", name + "_" + entry.second, params);
  }
}

void NEMLMechanicsAction::_add_scalar_aux(std::string name) {
  auto params = _factory.getValidParams("MaterialRealAux");

  params.set<MaterialPropertyName>("property") = name;
  params.set<AuxVariableName>("variable") = name;

  _problem->addAuxKernel("MaterialRealAux", name, params);
}
