#pragma once

#include "Action.h"

class NEMLMechanicsAction : public Action {
public:
  static InputParameters validParams();
  NEMLMechanicsAction(const InputParameters &params);

  virtual void act();

protected:
  void _add_tensor_variable(std::string name);
  void _add_scalar_variable(std::string name);

  void _add_tensor_aux(std::string name);
  void _add_scalar_aux(std::string name);

  std::vector<VariableName> _displacements;
  unsigned int _ndisp;
  bool _add_disp;
  bool _add_all;

  enum class Kinematics { Small, Large } _kinematics;
  enum class Formulation {Updated, Total} _formulation;

  std::map<Kinematics, bool> _kin_mapper = {{Kinematics::Small, false},
                                            {Kinematics::Large, true}};

  std::vector<MaterialPropertyName> _eigenstrains;
  std::vector<SubdomainName> _block;

  bool _homogenize;
  // Helper to translate into MOOSE talk
  const std::map<unsigned int, std::string> _order_mapper = 
  { {1, "FIRST"},
    {3, "THIRD"},
    {4, "FOURTH"},
    {6, "SIXTH"},
    {9, "NINTH"}};
  // Name of the homogenization scalar variable
  const std::string _hname = "hvar";
  // Name of the integrator
  const std::string _integrator_name = "integrator";
  // Other homogenization info
  std::vector<std::string> _constraint_types;
  std::vector<FunctionName> _targets;
};
