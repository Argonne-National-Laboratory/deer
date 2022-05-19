#include "InductionFunction.h"

registerMooseObject("DeerApp", InductionFunction);

InputParameters
InductionFunction::validParams()
{
  InputParameters params = Function::validParams();

  params.addRequiredParam<Point>("center", "Center point for calculating r");
  params.addRequiredParam<Point>("axis", "Tube axis");
  params.addRequiredParam<Real>("outer_radius", "Tube outer radius");
  params.addRequiredParam<Real>("thickness", "Tube thickness");
  params.addRequiredParam<Real>("delta", "Skin depth");

  return params;
};

InductionFunction::InductionFunction(const InputParameters & parameters)
  : Function(parameters),
    FunctionInterface(this),
    _center(getParam<Point>("center")),
    _axis(getParam<Point>("axis")),
    _radius(getParam<Real>("outer_radius")),
    _thickness(getParam<Real>("thickness")),
    _delta(getParam<Real>("delta"))
{
}

Real
InductionFunction::value(Real t, const Point & p) const
{
  Real r = (p - _center).cross(_axis).norm() / _axis.norm();
  Real x = _radius - r;

  return std::exp(-2 * x / _delta) /
         (M_PI * _radius * _delta * (1.0 - std::exp(-2 * _thickness / _delta)));
}
