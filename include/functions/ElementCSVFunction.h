#pragma once

#include "Function.h"

#include "DelimitedFileReader.h"

class ElementCSVFunction : public Function
{
public:
  static InputParameters validParams();
  ElementCSVFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const;

private:
  void read();

private:
  const std::string _file_name;
  MooseUtils::DelimitedFileReader _reader;
  const Real & _row;
  std::vector<std::vector<Real>> _data;
};
