#pragma once

#include "miconosstype.h"
#include <map>

class NLParameter {
public:
  NLParameter(const std::string &param_name, const double value);

  double getValue() const { return _value; }
  std::string getName() const { return _name; }
  void setValue(const double value) { _value = value; }
  NLParameter(const std::string &param_name, const double value,
              const bool rate_param);

  bool isRateParam() const { return _rate_parameter; };
  double getIncrement(const double dt) const;
  double getDRateDIncrement(const double dt) const;

protected:
  const std::string _name;
  double _value;
  const double _rate_parameter;
};

class NLSystemParameters {
public:
  NLSystemParameters(const std::vector<NLParameter *> params);
  NLSystemParameters(std::vector<std::string> param_names,
                     std::vector<double> param_values,
                     std::vector<std::string> rate_param_names,
                     vecD rate_param_values);

  uint getNParams() const { return _n_params; }
  uint getParamIndex(const std::string &pname) const;

  bool isRateParam(const uint index) const {
    return _params[index]->isRateParam();
  }
  bool isRateParam(const std::string &pname) const {
    return _params[getParamIndex(pname)]->isRateParam();
  }

  double getValue(const uint index) const { return _params[index]->getValue(); }
  double getValue(const std::string &pname) const {
    return _params[getParamIndex(pname)]->getValue();
  }

  double getIncrement(const uint index) const {
    return _params[index]->getIncrement(getValue("dt"));
  };
  double getIncrement(const std::string &pname) const {
    return _params[getParamIndex(pname)]->getIncrement(getValue("dt"));
  };

  double getDRateDIncrement(const uint index) const {
    return _params[index]->getIncrement(getDRateDIncrement("dt"));
  };
  double getDRateDIncrement(const std::string &pname) const {
    return _params[getParamIndex(pname)]->getDRateDIncrement(getValue("dt"));
  };

  void setValue(const uint index, const double &p) const {
    _params[index]->setValue(p);
  }
  void setValue(const std::string &pname, const double &p) const {
    _params[getParamIndex(pname)]->setValue(p);
  }

protected:
  std::vector<NLParameter *> _params;
  std::vector<NLParameter> _params_vector;
  std::map<std::string, uint> _name_index_map;
  const uint _n_params;
};
