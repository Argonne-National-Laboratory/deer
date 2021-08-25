#include "NLParameter.h"

NLParameter::NLParameter(const std::string &param_name, const double value)
    : _name(param_name), _value(value), _rate_parameter(false) {}

NLParameter::NLParameter(const std::string &param_name, const double value,
                         const bool rate_param)
    : _name(param_name), _value(value), _rate_parameter(rate_param) {}

double NLParameter::getValue() const { return _value; }
std::string NLParameter::getName() const { return _name; }
void NLParameter::setValue(const double value) { _value = value; }
bool NLParameter::isRateParam() const { return _rate_parameter; };
double NLParameter::getIncrement(const double dt) const {
  if (_rate_parameter)
    return getValue() * dt;
  else
    throw std::runtime_error(
        _name + " is not a rate parameter, you can't use getIncrement ");
};
double NLParameter::getDRateDIncrement(const double dt) const {
  if (_rate_parameter)
    return 1. / dt;
  else
    throw std::runtime_error(
        _name + " is not a rate parameter, you can't use getDRateDIncrement ");
};

NLSystemParameters::NLSystemParameters(std::vector<NLParameter *> params)
    : _params(params), _n_params(_params.size()) {
  for (uint i = 0; i < _n_params; i++)
    _name_index_map.insert(
        std::pair<std::string, uint>(_params[i]->getName(), i));
}

NLSystemParameters::NLSystemParameters(
    std::vector<std::string> param_names, vecD param_values,
    std::vector<std::string> rate_param_names, vecD rate_param_values)
    : _n_params(param_names.size() + rate_param_names.size()) {
  if (param_names.size() != param_values.size())
    throw std::runtime_error(
        "number of parameters and number of values do not match");

  if (rate_param_names.size() != rate_param_values.size())
    throw std::runtime_error(
        "number of rate parameters and number of values do not match");

  for (uint i = 0; i < param_names.size(); i++) {
    _name_index_map.insert(std::pair<std::string, uint>(param_names[i], i));
    _params_vector.push_back(NLParameter(param_names[i], param_values[i]));
  }

  uint k = 0;
  for (uint i = param_names.size(); i < _n_params; i++) {
    _name_index_map.insert(
        std::pair<std::string, uint>(rate_param_names[k], i));
    _params_vector.push_back(NLParameter(
        rate_param_names[k], rate_param_values[k], /*rate_param=*/true));
    k += 1;
  }

  for (uint i = 0; i < _n_params; i++)
    _params.push_back(&_params_vector[i]);
}

uint NLSystemParameters::getNParams() const { return _n_params; }

uint NLSystemParameters::getParamIndex(const std::string &pname) const {

  auto it = _name_index_map.find(pname);
  if (it == _name_index_map.end())
    throw std::runtime_error("can't find parameter " + pname +
                             " in the paramter list");

  return it->second;
}

bool NLSystemParameters::isRateParam(const uint index) const {
  return _params[index]->isRateParam();
}
bool NLSystemParameters::isRateParam(const std::string &pname) const {
  return _params[getParamIndex(pname)]->isRateParam();
}

double NLSystemParameters::getValue(const uint index) const {
  return _params[index]->getValue();
}
double NLSystemParameters::getValue(const std::string &pname) const {
  return _params[getParamIndex(pname)]->getValue();
}

double NLSystemParameters::getIncrement(const uint index) const {
  return _params[index]->getIncrement(getValue("dt"));
};
double NLSystemParameters::getIncrement(const std::string &pname) const {
  return _params[getParamIndex(pname)]->getIncrement(getValue("dt"));
};

double NLSystemParameters::getDRateDIncrement(const uint index) const {
  return _params[index]->getIncrement(getDRateDIncrement("dt"));
};
double NLSystemParameters::getDRateDIncrement(const std::string &pname) const {
  return _params[getParamIndex(pname)]->getDRateDIncrement(getValue("dt"));
};

void NLSystemParameters::setValue(const uint index, const double &p) const {
  _params[index]->setValue(p);
}
void NLSystemParameters::setValue(const std::string &pname,
                                  const double &p) const {
  _params[getParamIndex(pname)]->setValue(p);
}
