#include "NLVar.h"

NLVar::NLVar(const uint index, const std::string &var_name, const double x,
             const double x_old, const double scaling_factor)
    : _x(x), _x_old(x_old), _index(index), _scaling_factor(scaling_factor),
      _var_name(var_name) {}

double NLVar::realToNormalized(const double x) const {
  return x / _scaling_factor;
}

double NLVar::normalizedToReal(const double x) const {
  return x * _scaling_factor;
}

double NLVar::dRealdNormalized() const { return _scaling_factor; }
double NLVar::dNormalizeddReal() const { return 1. / _scaling_factor; }

double NLVar::getValue() const { return _x; }
double NLVar::getValueScaled() const { return realToNormalized(_x); }
double NLVar::getValueOld() const { return _x_old; }
double NLVar::getValueOldScaled() const { return realToNormalized(_x_old); }
double NLVar::getScalingFactor() const { return _scaling_factor; }
double NLVar::getDVarDVarScaled() const { return dRealdNormalized(); }
double NLVar::getDVarScaledDVar() const { return dNormalizeddReal(); }
double NLVar::getValueImplicit(const bool implicit) const {
  return implicit ? getValue() : getValueOld();
}

uint NLVar::getIndex() const { return _index; }
std::string NLVar::getName() const { return _var_name; }

void NLVar::setValue(const double &x) { _x = x; }
void NLVar::setValueFromScaled(const double &x) { _x = normalizedToReal(x); }
void NLVar::setValueOld(const double &x_old) { _x_old = x_old; }
void NLVar::setValueOldFromScaled(const double &x_old) {
  _x_old = normalizedToReal(_x_old);
}

void NLVar::setToOld() { setValue(_x_old); }

void NLVar::updateOldToCurrent() { setValueOld(_x); }

void NLVar::setValues(const double &x, const double &x_old) {
  setValue(x);
  setValueOld(x_old);
}

void NLVar::setScaleFactor(const double &sf) { _scaling_factor = sf; }

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

/*****************************************************************************
                                  GET METHODS
****************************************************************************/
double NLSystemVars::getValue(const uint &index) const {
  return _vars[index]->getValue();
}

double NLSystemVars::getValue(const std::string &vname) const {
  return _vars[getVarIndex(vname)]->getValue();
}

double NLSystemVars::getValueImplicit(const uint &index,
                                      const bool implicit) const {
  return _vars[index]->getValueImplicit(implicit);
}

double NLSystemVars::getValueImplicit(const std::string &vname,
                                      const bool implicit) const {
  return _vars[getVarIndex(vname)]->getValueImplicit(implicit);
}

double NLSystemVars::getValueScaled(const uint &index) const {
  return _vars[index]->getValueScaled();
}

double NLSystemVars::getValueScaled(const std::string &vname) const {
  return _vars[getVarIndex(vname)]->getValueScaled();
}

double NLSystemVars::getDVarScaledDVar(const uint &index) const {
  return _vars[index]->getDVarScaledDVar();
}

double NLSystemVars::getDVarScaledDVar(const std::string &vname) const {
  return _vars[getVarIndex(vname)]->getDVarScaledDVar();
}

double NLSystemVars::getDVarDVarScaled(const uint &index) const {
  return _vars[index]->getDVarDVarScaled();
}

double NLSystemVars::getDVarDVarScaled(const std::string &vname) const {
  return _vars[getVarIndex(vname)]->getDVarDVarScaled();
}

double NLSystemVars::getValueOld(const uint &index) const {
  return _vars[index]->getValueOld();
}

double NLSystemVars::getValueOld(const std::string &vname) const {
  return _vars[getVarIndex(vname)]->getValueOld();
}

double NLSystemVars::getValueOldScaled(const uint &index) const {
  return _vars[index]->getValueOldScaled();
}

double NLSystemVars::getValueOldScaled(const std::string &vname) const {
  return _vars[getVarIndex(vname)]->getValueOldScaled();
}

double NLSystemVars::getScalingFactor(const uint &index) const {
  return _vars[index]->getScalingFactor();
}

double NLSystemVars::getScalingFactor(const std::string &vname) const {
  return _vars[getVarIndex(vname)]->getScalingFactor();
}

std::string NLSystemVars::getName(const uint &index) const {
  return _vars[index]->getName();
}

/*****************************************************************************
                                  SET METHODS
****************************************************************************/

void NLSystemVars::setValue(const uint &index, const double &x) {
  _vars[index]->setValue(x);
}

void NLSystemVars::setValue(const std::string &vname, const double &x) {
  _vars[getVarIndex(vname)]->setValue(x);
}

void NLSystemVars::setValueFromScaled(const uint &index, const double &x) {
  _vars[index]->setValueFromScaled(x);
}

void NLSystemVars::setValueFromScaled(const std::string &vname,
                                      const double &x) {
  _vars[getVarIndex(vname)]->setValueFromScaled(x);
}

void NLSystemVars::setValueOld(const uint &index, const double &x_old) {
  _vars[index]->setValueOld(x_old);
}

void NLSystemVars::setValueOld(const std::string &vname, const double &x_old) {
  _vars[getVarIndex(vname)]->setValueOld(x_old);
}

void NLSystemVars::setValueOldFromScaled(const uint &index,
                                         const double &x_old) {
  _vars[index]->setValueOldFromScaled(x_old);
}

void NLSystemVars::setValueOldFromScaled(const std::string &vname,
                                         const double &x_old) {
  _vars[getVarIndex(vname)]->setValueOldFromScaled(x_old);
}

void NLSystemVars::setScaleFactor(const uint &index, const double &sf) {
  _vars[index]->setScaleFactor(sf);
}

void NLSystemVars::setScaleFactor(const std::string &vname, const double &sf) {
  _vars[getVarIndex(vname)]->setScaleFactor(sf);
}

void NLSystemVars::setToOld() {
  for (uint i = 0; i < _n_vars; i++)
    _vars[i]->setToOld();
}

void NLSystemVars::updateOldToCurrent() {
  for (uint i = 0; i < _n_vars; i++)
    _vars[i]->updateOldToCurrent();
}

void NLSystemVars::setFromVector(const vecD &new_value) {
  for (uint i = 0; i < _n_vars; i++)
    _vars[i]->setValue(new_value[i]);
}

void NLSystemVars::setOldFromVector(const vecD &new_value) {
  for (uint i = 0; i < _n_vars; i++)
    _vars[i]->setValueOld(new_value[i]);
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
