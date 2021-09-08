
#pragma once

#include "miconosstype.h"
#include <map>

/**definition of a miconoss non linear variable. A non linear variable is what
 * miconoss solves for. Each non linear variable requires an index, and a
 * name. The index identifies its position in the nonlinear system, the name is
 * used from getter and set methods of miconoss equations.
 **/
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

  /// set the variable current or old value
  ///{@
  void setValue(const double &x);
  void setValueOld(const double &x_old);
  ///@}

  /// set the variable current or old value given its scaled value
  ///{@
  void setValueFromScaled(const double &x);
  void setValueOldFromScaled(const double &x_old);
  ///@}

  /// set the current variable value to its old value
  void setToOld();
  /// set the old variable value to the current value
  void updateOldToCurrent();
  /// set both the current and old value
  void setValues(const double &x, const double &x_old);
  /// set the variable scaling factor
  void setScaleFactor(const double &sf);

protected:
  /// the current variable value
  double _x;
  /// the old variable value
  double _x_old;
  /// the variable index
  const uint _index;
  /// the variable scaling factor
  double _scaling_factor;
  /// the variable name
  const std::string _var_name;
};

/// a class assembling all non linear system used to solve a non linear system.
/// This class is used by NLSystem to solve the constrained or unconstrained non
/// linear problem.
class NLSystemVars {
public:
  /// constructor starting from variable pointers
  NLSystemVars(std::vector<NLVar *> vars);

  /// constructor starting from variables names, less flexible for intialization
  /// but take out the burden of allocating each single variable and providing
  /// consistent indeces. This method also hides variables pointer from the main
  /// program
  NLSystemVars(std::vector<std::string> var_names);

  /// get how many var the NL systems have
  uint getNVars() const { return _n_vars; }
  /// get teh var index given the var name
  uint getVarIndex(const std::string &vname) const;

  /*****************************************************************************
                                    GET METHODS
  ****************************************************************************/
  /// get a variable current value by index or name
  ///{@
  double getValue(const uint &index) const;
  double getValue(const std::string &vname) const;
  ///@}

  /// get a variable current implicit or explict value by index or name
  ///{@
  double getValueImplicit(const uint &index, const bool implicit) const;
  double getValueImplicit(const std::string &vname, const bool implicit) const;
  ///@}

  /// get a variable current scaled value by index or name
  ///{@
  double getValueScaled(const uint &index) const;
  double getValueScaled(const std::string &vname) const;
  ///@}

  /// get a variable dVar/dscaledVar by index or name
  ///{@
  double getDVarScaledDVar(const uint &index) const;
  double getDVarScaledDVar(const std::string &vname) const;
  ///@}

  /// get a variable dscaledVar/dVar by index or name
  ///{@
  double getDVarDVarScaled(const uint &index) const;
  double getDVarDVarScaled(const std::string &vname) const;
  ///@}

  /// get a variable old value by index or name
  ///{@
  double getValueOld(const uint &index) const;
  double getValueOld(const std::string &vname) const;
  ///@}

  /// get a variable scaled old value by index or name
  ///{@
  double getValueOldScaled(const uint &index) const;
  double getValueOldScaled(const std::string &vname) const;
  ///@}

  /// get a variable scaling factor by index or name
  ///{@
  double getScalingFactor(const uint &index) const;
  double getScalingFactor(const std::string &vname) const;
  ///@}

  /// get a variable name given the index
  std::string getName(const uint &index) const;

  /*****************************************************************************
                                    SET METHODS
  ****************************************************************************/

  /// set a variable current value by index or name
  ///{@
  void setValue(const uint &index, const double &x);
  void setValue(const std::string &vname, const double &x);
  ///@}

  /// set a variable value by index or name using a scaled value
  ///{@
  void setValueFromScaled(const uint &index, const double &x);
  void setValueFromScaled(const std::string &vname, const double &x);
  ///@}

  /// set a variable old value by index or name
  ///{@
  void setValueOld(const uint &index, const double &x_old);
  void setValueOld(const std::string &vname, const double &x_old);
  ///@}

  /// set a variable old value by index or name using a scaled value
  ///{@
  void setValueOldFromScaled(const uint &index, const double &x_old);
  void setValueOldFromScaled(const std::string &vname, const double &x_old);
  ///@}

  /// get a variable scaling factor by index or name
  ///{@
  void setScaleFactor(const uint &index, const double &sf);
  void setScaleFactor(const std::string &vname, const double &sf);
  ///@}

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
  /// the vector of NLVars pointer
  std::vector<NLVar *> _vars;
  /// the vector cof NLVars. This used when this class also initalize the
  /// variables
  std::vector<NLVar> _vars_vector;

  /// the name to index map
  std::map<std::string, uint> _name_index_map;
  /// the number of non linear variable
  const uint _n_vars;
};
