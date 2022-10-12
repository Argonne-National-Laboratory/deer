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
registerMooseAction("DeerApp", NEMLMechanicsAction, "add_scalar_kernel");
registerMooseAction("DeerApp", NEMLMechanicsAction, "add_user_object");

const std::vector<std::string> all_tensors = {
    "mechanical_strain", "stress", "elastic_strain", "inelastic_strain"};
const std::vector<std::string> all_scalars = {"energy", "dissipation"};

const std::vector<std::string> internal_strain_tensors = {"mechanical_strain_unrotated",
                                                          "inelastic_strain_unrotated"};

const std::map<std::pair<int, int>, std::string> tensor_map = {{std::make_pair(0, 0), "x-x"},
                                                               {std::make_pair(1, 1), "y-y"},
                                                               {std::make_pair(2, 2), "z-z"},
                                                               {std::make_pair(0, 1), "x-y"},
                                                               {std::make_pair(0, 2), "x-z"},
                                                               {std::make_pair(1, 2), "y-z"}};

InputParameters
NEMLMechanicsAction::validParams()
{
  InputParameters params = Action::validParams();

  params.addRequiredParam<std::vector<VariableName>>("displacements", "The displacement variables");
  params.addParam<bool>("add_displacements", true, "Add the displacement variables");

  MooseEnum kinematicType("small large", "small");
  params.addParam<MooseEnum>("kinematics", kinematicType, "Kinematic formulation");

  MooseEnum formulationType("updated total", "updated");
  params.addParam<MooseEnum>("formulation",
                             formulationType,
                             "Equilibrium formulation: updated or total "
                             "Lagrangian");

  params.addParam<bool>(
      "add_all_output", false, "Dump all the usual stress and strain variables to the output");

  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrains", std::vector<MaterialPropertyName>(), "Names of the eigenstrains");

  params.addParam<std::vector<SubdomainName>>(
      "block",
      "The list of subdomain names where neml mehcanisc should be used, "
      "default all blocks.");

  params.addParam<bool>("homogenize", false, "Apply homogenization constraints to the system.");
  params.addParam<std::vector<std::string>>("constraint_types",
                                            "Type of each constraint: "
                                            "stress or strain.");
  params.addParam<std::vector<FunctionName>>("targets",
                                             "Functions giving the target "
                                             "values of each constraint.");
  params.addParam<bool>("output_internal_strain",
                        false,
                        "add the internal inelastic and mechanical strain to the output");
  return params;
}

NEMLMechanicsAction::NEMLMechanicsAction(const InputParameters & params)
  : Action(params),
    _displacements(getParam<std::vector<VariableName>>("displacements")),
    _ndisp(_displacements.size()),
    _add_disp(getParam<bool>("add_displacements")),
    _add_all(getParam<bool>("add_all_output")),
    _kinematics(getParam<MooseEnum>("kinematics").getEnum<Kinematics>()),
    _formulation(getParam<MooseEnum>("formulation").getEnum<Formulation>()),
    _eigenstrains(getParam<std::vector<MaterialPropertyName>>("eigenstrains")),
    _block(params.isParamSetByUser("block") ? getParam<std::vector<SubdomainName>>("block")
                                            : std::vector<SubdomainName>(0)),
    _homogenize(getParam<bool>("homogenize")),
    _constraint_types(getParam<std::vector<std::string>>("constraint_types")),
    _targets(getParam<std::vector<FunctionName>>("targets"))
{
}

void
NEMLMechanicsAction::act()
{
  if (_current_task == "add_variable")
  {
    const bool second = _problem->mesh().hasSecondOrderElements();

    InputParameters params = _factory.getValidParams("MooseVariableBase");
    params.set<MooseEnum>("family") = "LAGRANGE";
    if (second)
      params.set<MooseEnum>("order") = "SECOND";
    else
      params.set<MooseEnum>("order") = "FIRST";

    for (const auto & disp_name : _displacements)
    {
      auto fe_type = AddVariableAction::feType(params);
      auto var_type = AddVariableAction::determineType(fe_type, 1);
      _problem->addVariable(var_type, disp_name, params);
    }

    if (_homogenize)
    {
      InputParameters params = _factory.getValidParams("MooseVariableBase");
      params.set<MooseEnum>("family") = "SCALAR";
      int order = 0;
      for (const auto & constraint_type : _constraint_types)
        if (MooseUtils::toUpper(constraint_type) != "NONE")
          order++;
      params.set<MooseEnum>("order") = order;
      auto fe_type = AddVariableAction::feType(params);
      auto var_type = AddVariableAction::determineType(fe_type, 1);
      _problem->addVariable(var_type, _hname, params);
    }
  }
  else if (_current_task == "add_material")
  {
    // Add the strain calculator
    auto params = _factory.getValidParams("ComputeNEMLStrain");

    params.set<std::vector<VariableName>>("displacements") = _displacements;
    params.set<std::vector<MaterialPropertyName>>("eigenstrain_names") = _eigenstrains;
    params.set<bool>("large_kinematics") = _kin_mapper[_kinematics];
    if (_homogenize)
      params.set<std::vector<VariableName>>("macro_gradient") = {_hname};

    _problem->addMaterial("ComputeNEMLStrain", "strain", params);
  }
  else if (_current_task == "add_kernel")
  {
    // Error check
    if ((_homogenize) && (_formulation == Formulation::Updated))
      mooseError("Homogenization constraints must be used with the total "
                 "lagrangian formulation");

    // Add the kernels
    for (unsigned int i = 0; i < _ndisp; ++i)
    {
      if (_formulation == Formulation::Updated)
      {
        auto params = _factory.getValidParams("StressDivergenceNEML");

        params.set<std::vector<VariableName>>("displacements") = _displacements;
        params.set<NonlinearVariableName>("variable") = _displacements[i];
        params.set<unsigned int>("component") = i;
        params.set<bool>("use_displaced_mesh") = _kin_mapper[_kinematics];
        if (_block.size() > 0)
          params.set<std::vector<SubdomainName>>("block") = _block;

        std::string name = "SD_" + Moose::stringify(i);

        _problem->addKernel("StressDivergenceNEML", name, params);
      }
      else if (_formulation == Formulation::Total)
      {
        auto params = _factory.getValidParams("TotalStressDivergenceNEML");

        params.set<std::vector<VariableName>>("displacements") = _displacements;
        params.set<NonlinearVariableName>("variable") = _displacements[i];
        params.set<unsigned int>("component") = i;
        params.set<bool>("large_kinematics") = _kin_mapper[_kinematics];
        if (_block.size() > 0)
          params.set<std::vector<SubdomainName>>("block") = _block;

        if (_homogenize)
        {
          params.set<std::vector<VariableName>>("macro_gradient") = {_hname};
          params.set<std::vector<std::string>>("constraint_types") = _constraint_types;
        }

        std::string name = "SD_" + Moose::stringify(i);

        _problem->addKernel("TotalStressDivergenceNEML", name, params);
      }
      else
      {
        mooseError("Unknown formulation type supplied to NEMLMechanics "
                   "action!");
      }
    }
  }
  else if (_current_task == "add_aux_variable")
  {
    if (_add_all)
    {
      for (auto name : all_tensors)
        _add_tensor_variable(name);
      if (getParam<bool>("output_internal_strain"))
        for (auto name : internal_strain_tensors)
          _add_tensor_variable(name);
      for (auto name : all_scalars)
        _add_scalar_variable(name);
    }
  }
  else if (_current_task == "add_aux_kernel")
  {
    if (_add_all)
    {
      for (auto name : all_tensors)
        _add_tensor_aux(name);
      if (getParam<bool>("output_internal_strain"))
        for (auto name : internal_strain_tensors)
          _add_tensor_aux(name);
      for (auto name : all_scalars)
        _add_scalar_aux(name);
    }
  }
  else if (_current_task == "add_scalar_kernel")
  {
    if (_homogenize)
    {
      InputParameters params = _factory.getValidParams("HomogenizationConstraintScalarKernel");
      params.set<NonlinearVariableName>("variable") = _hname;
      params.set<UserObjectName>("homogenization_constraint") = _integrator_name;
      if (_block.size() > 0)
        params.set<std::vector<SubdomainName>>("block") = _block;

      _problem->addScalarKernel(
          "HomogenizationConstraintScalarKernel", "HomogenizationConstraints", params);
    }
  }
  else if (_current_task == "add_user_object")
  {
    if (_homogenize)
    {
      InputParameters params = _factory.getValidParams("HomogenizationConstraint");
      params.set<std::vector<std::string>>("constraint_types") = _constraint_types;
      params.set<std::vector<FunctionName>>("targets") = _targets;
      params.set<bool>("large_kinematics") = _kin_mapper[_kinematics];
      params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_LINEAR};

      _problem->addUserObject("HomogenizationConstraint", _integrator_name, params);
    }
  }
}

void
NEMLMechanicsAction::_add_tensor_variable(std::string name)
{
  for (auto entry : tensor_map)
  {
    _add_scalar_variable(name + "_" + entry.second);
  }
}

void
NEMLMechanicsAction::_add_scalar_variable(std::string name)
{

  InputParameters params = _factory.getValidParams("MooseVariableBase");
  params.set<MooseEnum>("family") = "MONOMIAL";
  params.set<MooseEnum>("order") = "CONSTANT";

  auto fe_type = AddVariableAction::feType(params);
  auto var_type = AddVariableAction::determineType(fe_type, 1);

  _problem->addAuxVariable(var_type, name, params);
}

void
NEMLMechanicsAction::_add_tensor_aux(std::string name)
{
  for (auto entry : tensor_map)
  {
    auto params = _factory.getValidParams("RankTwoAux");

    params.set<MaterialPropertyName>("rank_two_tensor") = name;
    params.set<AuxVariableName>("variable") = name + "_" + entry.second;
    params.set<unsigned int>("index_i") = entry.first.first;
    params.set<unsigned int>("index_j") = entry.first.second;
    if (_block.size() > 0)
      params.set<std::vector<SubdomainName>>("block") = _block;

    _problem->addAuxKernel("RankTwoAux", name + "_" + entry.second, params);
  }
}

void
NEMLMechanicsAction::_add_scalar_aux(std::string name)
{
  auto params = _factory.getValidParams("MaterialRealAux");

  params.set<MaterialPropertyName>("property") = name;
  params.set<AuxVariableName>("variable") = name;
  if (_block.size() > 0)
    params.set<std::vector<SubdomainName>>("block") = _block;

  _problem->addAuxKernel("MaterialRealAux", name, params);
}
