#include "InterpolateCSVFunction.h"

registerMooseObject("DeerApp", InterpolateCSVFunction);

InputParameters
InterpolateCSVFunction::validParams()
{
  InputParameters params = Function::validParams();

  params.addRequiredParam<FileName>("filename", "CSV file name");
  params.addRequiredParam<Real>("row", "Row giving entry to use");
  params.declareControllable("row");

  params.addRequiredParam<Point>("a1", "First lattice vector");
  params.addParam<Point>("a2", Point(0, 1, 0), "Second lattice vector");
  params.addParam<Point>("a3", Point(0, 0, 1), "Third lattice vector");

  params.addRequiredParam<size_t>("n_x", "grid points in x");
  params.addParam<size_t>("n_y", 1, "grid points in y");
  params.addParam<size_t>("n_z", 1, "grad points in z");

  params.addParam<Real>("max_x", 1.0, "maximum x");
  params.addParam<Real>("max_y", 1.0, "maximum y");
  params.addParam<Real>("max_z", 1.0, "minimum z");

  params.addParam<Real>("min_x", 0.0, "minimum x");
  params.addParam<Real>("min_y", 0.0, "minimum y");
  params.addParam<Real>("min_z", 0.0, "minimum z");

  return params;
}

InterpolateCSVFunction::InterpolateCSVFunction(const InputParameters & parameters)
  : Function(parameters),
    _transform(RankTwoTensor(getParam<Point>("a1"), getParam<Point>("a2"), getParam<Point>("a3")).transpose().inverse()),
    _min(getParam<Real>("min_x"), getParam<Real>("min_y"), getParam<Real>("min_z")),
    _max(getParam<Real>("max_x"), getParam<Real>("max_y"), getParam<Real>("max_z")),
    _n({getParam<size_t>("n_x"), getParam<size_t>("n_y"), getParam<size_t>("n_z")}),
    _dim(_mci_feproblem.mesh().dimension()),
    _file_name(getParam<FileName>("filename")),
    _reader(_file_name),
    _row(getParam<Real>("row"))
{
  read();
}

void
InterpolateCSVFunction::read()
{
  _reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::ROWS);
  _reader.read();

  _data = _reader.getData();
}

Real
InterpolateCSVFunction::value(Real t, const Point & p) const
{
  // Transform back to original coordinates
  auto pp = _transform * p - _min;

  size_t id = 0;
  for (size_t i = 0; i < _dim; i++)
  {
    size_t ii = _dim - 1 - i;
    id *= _n[ii];
    Real h = (_max(ii) - _min(ii)) / _n[ii];
    size_t ind = floor(pp(ii) / h);
    id += ind;
  }

  size_t row = std::round(_row);

  if (row >= _data.size())
    mooseError("Requested row ", _row, " does not exist in CSV file data");

  if (_data[row].size() < id)
    mooseError("Row ", _row, " in CSV file does not have enough entries");

  return _data[row][id];
}
