#include "NLVar.h"

NLVar::NLVar(const uint index, const std::string &var_name, const double x,
             const double x_old, const double scaling_factor)
    : _x(x), _x_old(x_old), _index(index), _scaling_factor(scaling_factor),
      _var_name(var_name) {}

void NLVar::setValues(const double &x, const double &x_old) {
  setValue(x);
  setValueOld(x_old);
}

NLSystemVars::NLSystemVars(std::vector<NLVar *> vars)
    : _vars(vars), _n_vars(_vars.size()) {
  for (auto v : _vars)
    _name_index_map.insert(
        std::pair<std::string, uint>(v->getName(), v->getIndex()));
}

NLSystemVars::NLSystemVars(std::vector<std::string> var_names)
    : _n_vars(var_names.size()) {
  for (uint i = 0; i < _n_vars; i++) {
    _name_index_map.insert(std::pair<std::string, uint>(var_names[i], i));
    _vars_vector.push_back(NLVar(i, var_names[i]));
  }
  for (uint i = 0; i < _n_vars; i++)
    _vars.push_back(&_vars_vector[i]);
}

uint NLSystemVars::getVarIndex(const std::string &vname) const {

  auto it = _name_index_map.find(vname);
  if (it == _name_index_map.end())
    throw std::runtime_error("can't find variable " + vname + " in map");

  return it->second;
}

vecD NLSystemVars::getValueVector() const {
  vecD values(_n_vars);
  for (uint i = 0; i < _n_vars; i++)
    values[i] = getValue(i);
  return values;
}

vecD NLSystemVars::getValueVectorOld() const {
  vecD values(_n_vars);
  for (uint i = 0; i < _n_vars; i++)
    values[i] = getValueOld(i);
  return values;
}
