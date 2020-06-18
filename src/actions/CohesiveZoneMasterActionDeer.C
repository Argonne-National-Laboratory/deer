//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CohesiveZoneMasterActionDeer.h"
#include "Conversion.h"
#include "FEProblem.h"
#include "Factory.h"

registerMooseAction("DeerApp", CohesiveZoneMasterActionDeer,
                    "add_interface_kernel");

InputParameters CohesiveZoneMasterActionDeer::validParams() {
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Action to create an instance of the cohesive zone model kernel for "
      "each displacement component");
  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundary", "The list of boundary IDs from the mesh where the cohesive "
                  "zone will be applied");
  params.addRequiredParam<std::vector<VariableName>>(
      "displacements", "The displacements appropriate for the simulation "
                       "geometry and coordinate system");
  params.addParam<bool>("large_kinematics", true,
                        "if true(default) uses large kinematics formulation");
  return params;
}

CohesiveZoneMasterActionDeer::CohesiveZoneMasterActionDeer(
    const InputParameters &params)
    : Action(params) {}

void CohesiveZoneMasterActionDeer::act() {

  std::string kernel_name = getParam<bool>("large_kinematics")
                                ? "CZMInterfaceKernelPK1"
                                : "CZMInterfaceKernel";

  std::vector<VariableName> displacements;
  if (isParamValid("displacements"))
    displacements = getParam<std::vector<VariableName>>("displacements");

  if (_current_task == "add_interface_kernel") {
    for (unsigned int i = 0; i < displacements.size(); ++i) {
      // Create unique kernel name for each of the components
      std::string unique_kernel_name =
          kernel_name + "_" + _name + "_" + Moose::stringify(i);

      InputParameters paramsk = _factory.getValidParams(kernel_name);

      paramsk.set<bool>("use_displaced_mesh") = false;
      paramsk.set<unsigned int>("component") = i;
      paramsk.set<NonlinearVariableName>("variable") = displacements[i];
      paramsk.set<std::vector<VariableName>>("neighbor_var") = {
          displacements[i]};
      paramsk.set<std::vector<VariableName>>("displacements") = displacements;
      paramsk.set<std::vector<BoundaryName>>("boundary") =
          getParam<std::vector<BoundaryName>>("boundary");

      _problem->addInterfaceKernel(kernel_name, unique_kernel_name, paramsk);
    }
  }
}
