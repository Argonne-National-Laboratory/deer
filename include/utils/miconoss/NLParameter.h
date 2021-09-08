#pragma once

#include "miconosstype.h"
#include <map>

/*
This class implements a non linear system parameter, i.e. a quantity required to
fully define an equation. Paramters can also be provided in terms of their rate
 */
class NLParameter {
public:
  NLParameter(const std::string &param_name, const double value);
  NLParameter(const std::string &param_name, const double value,
              const bool rate_param);

  /// return a parameter value
  double getValue() const;
  /// return a parameter name
  std::string getName() const;
  /// set the parameter value
  void setValue(const double value);
  /// if true means the parameter is provided as a rate quantitiy
  bool isRateParam() const;
  /// return the increment of a rate parameter
  double getIncrement(const double dt) const;
  /// return the dpdot/dDeltaP
  double getDRateDIncrement(const double dt) const;

protected:
  const std::string _name;
  double _value;
  const double _rate_parameter;
};

/*
This class is a collection of parameters and is used by other classes toa ccess
all the requried parameters. It also provides methods to get and set parameters.
 */
class NLSystemParameters {
public:
  NLSystemParameters(const std::vector<NLParameter *> params);
  NLSystemParameters(std::vector<std::string> param_names,
                     std::vector<double> param_values,
                     std::vector<std::string> rate_param_names,
                     vecD rate_param_values);

  /// return the number of parameters
  uint getNParams() const;
  /// return the index of a specifc parameter given its name
  uint getParamIndex(const std::string &pname) const;

  /// return true if the parameter assocaite to index index is a rate parameter
  bool isRateParam(const uint index) const;
  /// return true if the parameter named pname is a rate parameter
  bool isRateParam(const std::string &pname) const;

  /// get a parameter value by index or name
  ///@{
  double getValue(const uint index) const;
  double getValue(const std::string &pname) const;
  ///@}

  /// get the increment of a rate parameter by index or name
  ///@{
  double getIncrement(const uint index) const;
  double getIncrement(const std::string &pname) const;
  ///@}

  /// get the derivative of a rate parameter w.r.t. its increment
  ///@{
  double getDRateDIncrement(const uint index) const;
  double getDRateDIncrement(const std::string &pname) const;
  ///@}

  /// set a parameter value by index or name
  ///@{
  void setValue(const uint index, const double &p) const;
  void setValue(const std::string &pname, const double &p) const;
  ///@}

protected:
  ///  vector of parameter pointers
  std::vector<NLParameter *> _params;
  ///  vector parameters, reuired if parameters are defined by this class
  std::vector<NLParameter> _params_vector;
  ///  map between paramter name and index
  std::map<std::string, uint> _name_index_map;
  /// the total number of parameters
  const uint _n_params;
};
