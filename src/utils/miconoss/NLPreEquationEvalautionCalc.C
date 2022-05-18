#include "NLPreEquationEvalautionCalc.h"

NLPreEquationEvalautionCalc::NLPreEquationEvalautionCalc(
    NLSystemVars * const sysvars,
    const NLSystemParameters * sysparams,
    const std::vector<std::string> & value_names)
  : _sys_vars(sysvars),
    _sysparams(sysparams),
    _n_values(value_names.size()),
    _n_vars(_sys_vars->getNVars()),
    _values(_n_values),
    _dvalues_dx(_n_values, std::vector<double>(_n_vars)),
    _values_explicit(_n_values),
    _dvalues_dx_explicit(_n_vars, 0)
{

  for (uint i = 0; i < _n_values; i++)
    _name_index_map.insert(std::pair<std::string, uint>(value_names[i], i));
}

void
NLPreEquationEvalautionCalc::updateAll(const bool updade_explicit)
{
  if (updade_explicit)
    updateValuesExplicit();

  updateValues();
  updateDerivatives();
}

uint
NLPreEquationEvalautionCalc::getValueIndex(const std::string & vname) const
{

  auto it = _name_index_map.find(vname);
  if (it == _name_index_map.end())
    throw std::runtime_error("can't find precalaulcate value " + vname +
                             " in the precalcualted parameter list");

  return it->second;
}

double
NLPreEquationEvalautionCalc::getValue(const std::string & vname, const bool implicit) const
{
  if (implicit)
    return _values[getValueIndex(vname)];
  else
    return _values_explicit[getValueIndex(vname)];
}

vecD
NLPreEquationEvalautionCalc::getDValueDX(const std::string & vname, const bool implicit) const
{
  if (implicit)
    return _dvalues_dx[getValueIndex(vname)];
  else
    return _dvalues_dx_explicit;
}

void
NLPreEquationEvalautionCalc::setValue(const std::string & vname, double value, const bool implicit)
{
  if (implicit)
    _values[getValueIndex(vname)] = value;
  else
    _values_explicit[getValueIndex(vname)] = value;
}

void
NLPreEquationEvalautionCalc::setDValueDX(const std::string & vname,
                                         const vecD & dvaluedx,
                                         const bool implicit)
{
  if (implicit)
    _dvalues_dx[getValueIndex(vname)] = dvaluedx;
}
