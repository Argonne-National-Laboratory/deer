#pragma once

#include "Function.h"

#include "DelimitedFileReader.h"

class InterpolateCSVFunction : public Function
{
public:
  static InputParameters validParams();
  InterpolateCSVFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const;

private:
  void read();

private:
  /// the three lattice vectors
  RankTwoTensor _transform;

  Point _min, _max;
  std::vector<size_t> _n;
  size_t _dim;

  const std::string _file_name;
  MooseUtils::DelimitedFileReader _reader;
  const Real & _row;
  std::vector<std::vector<Real>> _data;
};
