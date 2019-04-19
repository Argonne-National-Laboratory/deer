#include "NEMLMechanicsAction.h"

#include "FEProblem.h"
#include "MooseMesh.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("DeerApp", NEMLMechanicsAction, "add_variable");
registerMooseAction("DeerApp", NEMLMechanicsAction, "add_kernel");
registerMooseAction("DeerApp", NEMLMechanicsAction, "add_material");
registerMooseAction("DeerApp", NEMLMechanicsAction, "add_aux_variable");
registerMooseAction("DeerApp", NEMLMechanicsAction, "add_aux_kernel");

const std::vector<std::string> all_tensors = {"mechanical_strain", "stress",
  "elastic_strain", "inelastic_strain"};
const std::vector<std::string> all_scalars = {"energy", "dissipation"};

const std::map<std::pair<int,int>, std::string> tensor_map = {
  {std::make_pair(0,0), "x-x"},
  {std::make_pair(1,1), "y-y"},
  {std::make_pair(2,2), "z-z"},
  {std::make_pair(0,1), "x-y"},
  {std::make_pair(0,2), "x-z"},
  {std::make_pair(1,2), "y-z"}};

template <>
InputParameters
validParams<NEMLMechanicsAction>()
{
  InputParameters params = validParams<Action>();

  params.addRequiredParam<std::vector<VariableName>>(
      "displacements", "The displacement variables");
  params.addParam<bool>("add_displacements", true, "Add the displacement variables");

  MooseEnum kinematicType("small large", "small");
  params.addParam<MooseEnum>("kinematics", kinematicType, "Kinematic formulation");

  params.addParam<bool>("add_all_output", false, "Dump all the usual stress and strain variables to the output");

  return params;
}

NEMLMechanicsAction::NEMLMechanicsAction(const InputParameters & params) : 
    Action(params),
    _displacements(getParam<std::vector<VariableName>>("displacements")),
    _ndisp(_displacements.size()),
    _add_disp(getParam<bool>("add_displacements")),
    _add_all(getParam<bool>("add_all_output")),
    _kinematics(getParam<MooseEnum>("kinematics").getEnum<Kinematics>())
{

}

void 
NEMLMechanicsAction::act()
{
  if (_current_task == "add_variable")
  {
    const bool second = _problem->mesh().hasSecondOrderElements();

    for (const auto & disp : _displacements)
    {
      _problem->addVariable(disp,
                           FEType(Utility::string_to_enum<Order>(second ? "SECOND" : "FIRST"),
                                  Utility::string_to_enum<FEFamily>("LAGRANGE")),
                           1.0,
                           nullptr);
    }
  }
  else if (_current_task == "add_material")
  {
    // Add the strain calculator
    auto params = _factory.getValidParams("ComputeNEMLStrain");

    params.set<std::vector<VariableName>>("displacements") = _displacements;
    params.set<bool>("use_displaced_mesh") = _kin_mapper[_kinematics];

    _problem->addMaterial("ComputeNEMLStrain", "strain", params);
  }
  else if (_current_task == "add_kernel")
  {
    // Add the kernels
    for (unsigned int i = 0; i < _ndisp; ++i)
    {
      auto params = _factory.getValidParams("StressDivergenceNEML");

      params.set<std::vector<VariableName>>("displacements") = _displacements;
      params.set<NonlinearVariableName>("variable") = _displacements[i];
      params.set<unsigned int>("component") = i;
      params.set<bool>("use_displaced_mesh") = _kin_mapper[_kinematics];
      
      std::string name = "SD_" + Moose::stringify(i);

      _problem->addKernel("StressDivergenceNEML", name, params);
    }
  }
  else if (_current_task == "add_aux_variable") {
    if (_add_all) {
      for (auto name : all_tensors) _add_tensor_variable(name);
      for (auto name : all_scalars) _add_scalar_variable(name);
    }
  }
  else if (_current_task == "add_aux_kernel")
  {
    if (_add_all) {
      for (auto name : all_tensors) _add_tensor_aux(name);
      for (auto name : all_scalars) _add_scalar_aux(name);
    }
  }
}

void
NEMLMechanicsAction::_add_tensor_variable(std::string name)
{
  for (auto entry : tensor_map) {
    _add_scalar_variable(name + "_" + entry.second);
  }
}

void
NEMLMechanicsAction::_add_scalar_variable(std::string name)
{
  _problem->addAuxVariable(name, 
                           FEType(Utility::string_to_enum<Order>("CONSTANT"),
                                  Utility::string_to_enum<FEFamily>("MONOMIAL")),
                           nullptr);
}

void
NEMLMechanicsAction::_add_tensor_aux(std::string name)
{
  for (auto entry : tensor_map) {
    auto params = _factory.getValidParams("RankTwoAux");

    params.set<MaterialPropertyName>("rank_two_tensor") = name;
    params.set<AuxVariableName>("variable") = name + "_" + entry.second;
    params.set<unsigned int>("index_i") = entry.first.first;
    params.set<unsigned int>("index_j") = entry.first.second;

    _problem->addAuxKernel("RankTwoAux", name + "_" + entry.second, params);
  }
}

void
NEMLMechanicsAction::_add_scalar_aux(std::string name)
{
  auto params = _factory.getValidParams("MaterialRealAux");

  params.set<MaterialPropertyName>("property") = name;
  params.set<AuxVariableName>("variable") = name;

  _problem->addAuxKernel("MaterialRealAux", name, params);
}
