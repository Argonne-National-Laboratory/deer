
#pragma once

#include "miconosstype.h"
#include <map>

class NLVar {
public:
  NLVar(const uint index, const std::string &var_name, const double x = 0,
        const double x_old = 0, const double scaling_factor = 1);

  double getValue() const { return _x; }
  double getValueScaled() const { return _x / _scaling_factor; }
  double getValueOld() const { return _x_old; }
  double getValueOldScaled() const { return _x_old / _scaling_factor; }
  double getScalingFactor() const { return _scaling_factor; }
  double getDVarDVarScaled() const { return _scaling_factor; }
  double getDVarScaledDVar() const { return 1. / _scaling_factor; }
  uint getIndex() const { return _index; }
  std::string getName() const { return _var_name; }

  void setValue(const double &x) { _x = x; }
  void setValueFromScaled(const double &x) { _x = x * _scaling_factor; }
  void setValueOld(const double &x_old) { _x_old = x_old; }
  void setValueOldFromScaled(const double &x_old) {
    _x_old = x_old * _scaling_factor;
  }

  void setToOld() { setValue(_x_old); };
  void updateOldToCurrent() { setValueOld(_x); };
  void setValues(const double &x, const double &x_old);
  void setScaleFactor(const double &sf) { _scaling_factor = sf; }

  double getValueImplicit(const bool implicit) const {
    return implicit ? getValue() : getValueOld();
  }

protected:
  double _x;
  double _x_old;
  const uint _index;
  double _scaling_factor;
  const std::string _var_name;
};

class NLSystemVars {
public:
  /// constructor starting from varaible pointers
  NLSystemVars(std::vector<NLVar *> vars);

  /// constructor starting from variables names, less flexible for intialization
  /// but take out the burdne of allocating each single variable and providing
  /// consistne indeces. THis method also hides varaible from the main program
  NLSystemVars(std::vector<std::string> var_names);

  uint getNVars() const { return _n_vars; }
  uint getVarIndex(const std::string &vname) const;

  /*****************************************************************************
                                    GET METHODS
  ****************************************************************************/
  double getValue(const uint &index) const { return _vars[index]->getValue(); }
  double getValue(const std::string &vname) const {
    return _vars[getVarIndex(vname)]->getValue();
  }

  double getValueImplicit(const uint &index, const bool implicit) const {
    return _vars[index]->getValueImplicit(implicit);
  }
  double getValueImplicit(const std::string &vname, const bool implicit) const {
    return _vars[getVarIndex(vname)]->getValueImplicit(implicit);
  }

  double getValueScaled(const uint &index) const {
    return _vars[index]->getValueScaled();
  }
  double getValueScaled(const std::string &vname) const {
    return _vars[getVarIndex(vname)]->getValueScaled();
  }

  double getDVarScaledDVar(const uint &index) const {
    return _vars[index]->getDVarScaledDVar();
  }
  double getDVarScaledDVar(const std::string &vname) const {
    return _vars[getVarIndex(vname)]->getDVarScaledDVar();
  }

  double getDVarDVarScaled(const uint &index) const {
    return _vars[index]->getDVarDVarScaled();
  }
  double getDVarDVarScaled(const std::string &vname) const {
    return _vars[getVarIndex(vname)]->getDVarDVarScaled();
  }

  double getValueOld(const uint &index) const {
    return _vars[index]->getValueOld();
  }
  double getValueOld(const std::string &vname) const {
    return _vars[getVarIndex(vname)]->getValueOld();
  }
  double getValueOldScaled(const uint &index) const {
    return _vars[index]->getValueOldScaled();
  }
  double getValueOldScaled(const std::string &vname) const {
    return _vars[getVarIndex(vname)]->getValueOldScaled();
  }

  double getScalingFactor(const uint &index) const {
    return _vars[index]->getScalingFactor();
  }
  double getScalingFactor(const std::string &vname) const {
    return _vars[getVarIndex(vname)]->getScalingFactor();
  }

  std::string getName(const uint &index) const {
    return _vars[index]->getName();
  }

  /*****************************************************************************
                                    SET METHODS
  ****************************************************************************/
  void setValue(const uint &index, const double &x) {
    _vars[index]->setValue(x);
  }
  void setValue(const std::string &vname, const double &x) {
    _vars[getVarIndex(vname)]->setValue(x);
  }
  void setValueFromScaled(const uint &index, const double &x) {
    _vars[index]->setValueFromScaled(x);
  }
  void setValueFromScaled(const std::string &vname, const double &x) {
    _vars[getVarIndex(vname)]->setValueFromScaled(x);
  }

  void setValueOld(const uint &index, const double &x_old) {
    _vars[index]->setValueOld(x_old);
  }
  void setValueOld(const std::string &vname, const double &x_old) {
    _vars[getVarIndex(vname)]->setValueOld(x_old);
  }
  void setValueOldFromScaled(const uint &index, const double &x_old) {
    _vars[index]->setValueOldFromScaled(x_old);
  }
  void setValueOldFromScaled(const std::string &vname, const double &x_old) {
    _vars[getVarIndex(vname)]->setValueOldFromScaled(x_old);
  }

  void setScaleFactor(const uint &index, const double &sf) {
    _vars[index]->setScaleFactor(sf);
  }
  void setScaleFactor(const std::string &vname, const double &sf) {
    _vars[getVarIndex(vname)]->setScaleFactor(sf);
  }

  /**********************************************************************************
  update and reset methods for all variables, typically used after solving the
  the non linear sytem or while substepping.
  **********************************************************************************/

  void setToOld() {
    for (uint i = 0; i < _n_vars; i++)
      _vars[i]->setToOld();
  }
  void updateOldToCurrent() {
    for (uint i = 0; i < _n_vars; i++)
      _vars[i]->updateOldToCurrent();
  }
  void setFromVector(const vecD &new_value) {
    for (uint i = 0; i < _n_vars; i++)
      _vars[i]->setValue(new_value[i]);
  }
  void setOldFromVector(const vecD &new_value) {
    for (uint i = 0; i < _n_vars; i++)
      _vars[i]->setValueOld(new_value[i]);
  }

  /**********************************************************************************
  varaiable value maps, useful for rate equations where we need both values for
  using the theta intgeration method
  **********************************************************************************/
  vecD getValueVector() const;
  vecD getValueVectorOld() const;

protected:
  std::vector<NLVar *> _vars;
  std::vector<NLVar> _vars_vector;
  std::map<std::string, uint> _name_index_map;
  const uint _n_vars;
};
