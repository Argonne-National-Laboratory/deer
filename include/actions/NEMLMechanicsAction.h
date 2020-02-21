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

  std::map<Kinematics, bool> _kin_mapper = {{Kinematics::Small, false},
                                            {Kinematics::Large, true}};

  std::vector<MaterialPropertyName> _eigenstrains;
  std::vector<SubdomainName> _block;
};
