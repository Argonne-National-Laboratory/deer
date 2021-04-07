#pragma once

#include "miconosstype.h"
#include <map>

class NLParameter {
public:
  NLParameter(const std::string &param_name, const double value);
  NLParameter(const std::string &param_name, const double value,
              const bool rate_param);

  double getValue() const;
  std::string getName() const;
  void setValue(const double value);
  bool isRateParam() const;
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

  uint getNParams() const;
  uint getParamIndex(const std::string &pname) const;

  bool isRateParam(const uint index) const;
  bool isRateParam(const std::string &pname) const;

  double getValue(const uint index) const;
  double getValue(const std::string &pname) const;

  double getIncrement(const uint index) const;
  double getIncrement(const std::string &pname) const;

  double getDRateDIncrement(const uint index) const;
  double getDRateDIncrement(const std::string &pname) const;

  void setValue(const uint index, const double &p) const;
  void setValue(const std::string &pname, const double &p) const;

protected:
  std::vector<NLParameter *> _params;
  std::vector<NLParameter> _params_vector;
  std::map<std::string, uint> _name_index_map;
  const uint _n_params;
};
