
#pragma once

#include "miconosstype.h"
#include <map>

class NLVar {
public:
  NLVar(const uint index, const std::string &var_name, const double x = 0,
        const double x_old = 0, const double scaling_factor = 1);

  /// methods to override to achive custom scaling
  ///{@
  virtual double realToNormalized(const double x) const;
  virtual double normalizedToReal(const double x) const;
  virtual double dRealdNormalized() const;
  virtual double dNormalizeddReal() const;
  ///@}

  /// getter methods
  ///{@
  double getValue() const;
  double getValueScaled() const;
  double getValueOld() const;
  double getValueOldScaled() const;
  double getScalingFactor() const;
  double getDVarDVarScaled() const;
  double getDVarScaledDVar() const;
  double getValueImplicit(const bool implicit) const;
  uint getIndex() const;
  std::string getName() const;

  /// set methods
  ///{@
  void setValue(const double &x);
  void setValueFromScaled(const double &x);
  void setValueOld(const double &x_old);
  void setValueOldFromScaled(const double &x_old);
  ///@}

  void setToOld();
  void updateOldToCurrent();
  void setValues(const double &x, const double &x_old);
  void setScaleFactor(const double &sf);

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
  /// but take out the burden of allocating each single variable and providing
  /// consistent indeces. This method also hides variables pointer from the main
  /// program
  NLSystemVars(std::vector<std::string> var_names);

  uint getNVars() const { return _n_vars; }
  uint getVarIndex(const std::string &vname) const;

  /*****************************************************************************
                                    GET METHODS
  ****************************************************************************/
  double getValue(const uint &index) const;
  double getValue(const std::string &vname) const;
  double getValueImplicit(const uint &index, const bool implicit) const;
  double getValueImplicit(const std::string &vname, const bool implicit) const;
  double getValueScaled(const uint &index) const;
  double getValueScaled(const std::string &vname) const;
  double getDVarScaledDVar(const uint &index) const;
  double getDVarScaledDVar(const std::string &vname) const;
  double getDVarDVarScaled(const uint &index) const;
  double getDVarDVarScaled(const std::string &vname) const;
  double getValueOld(const uint &index) const;
  double getValueOld(const std::string &vname) const;
  double getValueOldScaled(const uint &index) const;
  double getValueOldScaled(const std::string &vname) const;
  double getScalingFactor(const uint &index) const;
  double getScalingFactor(const std::string &vname) const;
  std::string getName(const uint &index) const;

  /*****************************************************************************
                                    SET METHODS
  ****************************************************************************/
  void setValue(const uint &index, const double &x);
  void setValue(const std::string &vname, const double &x);
  void setValueFromScaled(const uint &index, const double &x);
  void setValueFromScaled(const std::string &vname, const double &x);
  void setValueOld(const uint &index, const double &x_old);
  void setValueOld(const std::string &vname, const double &x_old);
  void setValueOldFromScaled(const uint &index, const double &x_old);
  void setValueOldFromScaled(const std::string &vname, const double &x_old);
  void setScaleFactor(const uint &index, const double &sf);
  void setScaleFactor(const std::string &vname, const double &sf);

  /**********************************************************************************
  update and reset and get methods for all variables, typically used after
  solving the the non linear sytem or while substepping.
  **********************************************************************************/

  void setToOld();
  void updateOldToCurrent();
  void setFromVector(const vecD &new_value);
  void setOldFromVector(const vecD &new_value);
  vecD getValueVector() const;
  vecD getValueVectorOld() const;

protected:
  std::vector<NLVar *> _vars;
  std::vector<NLVar> _vars_vector;
  std::map<std::string, uint> _name_index_map;
  const uint _n_vars;
};
