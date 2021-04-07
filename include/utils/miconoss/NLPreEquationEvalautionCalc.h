#pragma once

#include "NLParameter.h"
#include "NLVar.h"
#include <map>

class NLPreEquationEvalautionCalc {
public:
  NLPreEquationEvalautionCalc(NLSystemVars *const sysvars,
                              const NLSystemParameters *sysparams,
                              const std::vector<std::string> &value_names);

  /**
   * methods to override filling values an their derivatives
   */
  /**@{*/
  virtual void updateValues(const bool implicit = true){};
  virtual void updateDerivatives(const bool implicit = true){};
  /**@}*/

  void updateValuesExplicit() { updateValues(/*implicit =*/false); };

  void updateAll(const bool updade_explicit = false);

  double getValue(const std::string &vname, const bool implicit = true) const;

  vecD getDValueDX(const std::string &vname, const bool implicit = true) const;

protected:
  uint getValueIndex(const std::string &vname) const;

  void setValue(const std::string &vname, double value,
                const bool implicit = true);

  void setDValueDX(const std::string &vname, const vecD &dvaluedx,
                   const bool implicit = true);

  NLSystemVars *const _sys_vars;
  const NLSystemParameters *_sysparams;
  const uint _n_values;
  const uint _n_vars;
  vecD _values;
  std::vector<vecD> _dvalues_dx;
  vecD _values_explicit;
  const vecD _dvalues_dx_explicit;
  std::map<std::string, uint> _name_index_map;
};
