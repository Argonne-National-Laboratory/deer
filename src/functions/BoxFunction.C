#include "BoxFunction.h"

registerMooseObject("DeerApp", BoxFunction);

InputParameters
BoxFunction::validParams()
{
  InputParameters params = Function::validParams();

  params.addParam<Point>("lower_bound", "One point of box");
  params.addParam<Point>("upper_bound", "Second point of box");
  params.addParam<Real>("value", 0.0, "Value inside box");

  return params;
}

BoxFunction::BoxFunction(const InputParameters & parameters)
  : Function(parameters),
    _lb(getParam<Point>("lower_bound")),
    _ub(getParam<Point>("upper_bound")),
    _value(getParam<Real>("value"))
{
}

Real
BoxFunction::value(Real t, const Point & p) const
{
  for (unsigned int i = 0; i < LIBMESH_DIM; i++)
    if ((p(i) < std::min(_lb(i), _ub(i))) || (p(i) > std::max(_ub(i), _lb(i))))
      return 0.0;

  return _value;
}
