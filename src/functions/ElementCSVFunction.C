#include "ElementCSVFunction.h"

registerMooseObject("DeerApp", ElementCSVFunction);

InputParameters
ElementCSVFunction::validParams()
{
  InputParameters params = Function::validParams();

  params.addRequiredParam<FileName>("filename", "CSV file name");
  params.addRequiredParam<Real>("row", "Row giving entry to use");
  params.declareControllable("row");

  return params;
}

ElementCSVFunction::ElementCSVFunction(const InputParameters & parameters)
  : Function(parameters),
    _file_name(getParam<FileName>("filename")),
    _reader(_file_name),
    _row(getParam<Real>("row"))
{
  read();
}

void
ElementCSVFunction::read()
{
  _reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::ROWS);
  _reader.read();
  
  _data = _reader.getData();
}

Real
ElementCSVFunction::value(Real t, const Point & p) const
{
  std::set<const Elem *> potential_elements;
  (*(_ti_feproblem.mesh().getPointLocator()))(p, potential_elements);
  if (potential_elements.size() == 0)
    mooseError("No element located at point ", p);
  else if (potential_elements.size() > 1)
    mooseWarning("More than one element located at point ", p, 
                 ".  Will use lowest element id");

  const Elem * use = nullptr;
  for (const auto & elem : potential_elements)
    if (!use || elem->id() < use->id())
      use = elem;
  
  size_t row = std::round(_row);

  if (row >= _data.size())
    mooseError("Requested row ", _row, " does not exist in CSV file data");

  if (_data[row].size() < _ti_feproblem.mesh().nElem())
    mooseError("Row ", _row, " in CSV file does not have enough entries");
  
  return _data[row][use->id()];
}
