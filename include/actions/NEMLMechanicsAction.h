#ifndef NEMLMECHANICSACTION_H
#define NEMLMECHANICSACTION_H

#include "Action.h"

class NEMLMechanicsAction;

template <>
InputParameters validParams<NEMLMechanicsAction>();

class NEMLMechanicsAction : public Action
{
 public:
  NEMLMechanicsAction(const InputParameters & params);

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

  enum class Kinematics
  {
    Small,
    Large
  } _kinematics;

  std::map<Kinematics, bool> _kin_mapper = {
    {Kinematics::Small, false},
    {Kinematics::Large, true}};
};

#endif // NEMLMECHANICSACTION_H
