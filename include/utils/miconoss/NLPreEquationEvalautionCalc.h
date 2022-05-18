#pragma once

#include "NLParameter.h"
#include "NLVar.h"
#include <map>

/*
Class defining an object producing variable dependent quantites that needs to,
or can be calcaulted precalcualted before updating the non linear system
equations. This is handy when some value is needed by more than one equation.
*/
class NLPreEquationEvalautionCalc
{
public:
  NLPreEquationEvalautionCalc(NLSystemVars * const sysvars,
                              const NLSystemParameters * sysparams,
                              const std::vector<std::string> & value_names);

  /**
   * methods to override filling values an their derivatives
   */

  /// method computing all the values
  virtual void updateValues(const bool implicit = true){};
  /// method computing the derivates of all the values
  virtual void updateDerivatives(const bool implicit = true){};

  /**
   * other methods of the class
   */
  /// update all values using an explicit update
  void updateValuesExplicit() { updateValues(/*implicit =*/false); };

  /// method updating all values
  void updateAll(const bool updade_explicit = false);

  /**
   * methods typically used by the equations to retrive precalculated values
   */
  /// return a precomputede value by name
  double getValue(const std::string & vname, const bool implicit = true) const;

  /// return a precomputede dvalue/dX by name
  vecD getDValueDX(const std::string & vname, const bool implicit = true) const;

protected:
  /// return the index of a value given its name
  uint getValueIndex(const std::string & vname) const;

  /**
   * methods typically used within updateValues and updateDerivatives
   */
  /// set a value given its name
  void setValue(const std::string & vname, double value, const bool implicit = true);
  /// set a value gradient given its name
  void setDValueDX(const std::string & vname, const vecD & dvaluedx, const bool implicit = true);

  /**
   * internal object required to store values
   */
  /// the nonlinear system variables
  NLSystemVars * const _sys_vars;
  /// the nonlinear system paramters
  const NLSystemParameters * _sysparams;
  /// the number of values computed by this object
  const uint _n_values;
  /// the number of variables of the non linear system
  const uint _n_vars;
  /// the vector of computed values
  vecD _values;
  /// the vector of values gradient
  std::vector<vecD> _dvalues_dx;
  /// the value calcualted using an explicit update
  vecD _values_explicit;
  /// the derivative of the explicit value
  const vecD _dvalues_dx_explicit;
  /// the map linking a value name to its index
  std::map<std::string, uint> _name_index_map;
};
