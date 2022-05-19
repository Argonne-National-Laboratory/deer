#include "OneOffBendingFunction.h"

registerMooseObject("DeerApp", OneOffBendingFunction);

InputParameters
OneOffBendingFunction::validParams()
{
  InputParameters params = Function::validParams();

  params.addRequiredParam<unsigned int>("component", "Which component: x (0), y (1), or z (2)");

  params.addRequiredParam<Real>("length", "Length of pipe");
  params.addRequiredParam<Real>("radius", "Radius of bend");
  params.addRequiredParam<FunctionName>("theta", "Theta as a function of time");

  return params;
}

OneOffBendingFunction::OneOffBendingFunction(const InputParameters & parameters)
  : Function(parameters),
    FunctionInterface(this),
    _component(getParam<unsigned int>("component")),
    _length(getParam<Real>("length")),
    _radius(getParam<Real>("radius")),
    _theta(getFunction("theta"))
{
}

Real
OneOffBendingFunction::value(Real t, const Point & pt) const
{
  return displacements(t, pt)(_component);
}

RealVectorValue
OneOffBendingFunction::displacements(Real t, const Point & pt) const
{
  RealVectorValue disp;
  Real theta = _theta.value(t, pt);

  Real L = _length - pt(2);
  Real dz = _radius * theta;

  Real h = pt(0);

  if (L > _radius * theta)
  {
    disp(0) = 0.0;
    disp(1) = 0.0;
    disp(2) = dz;
  }
  else
  {
    Real theta_p = theta - L / _radius;
    disp(0) = -(1.0 - cos(theta_p)) * _radius - h * (1 - cos(theta_p));
    disp(1) = 0.0;
    disp(2) = dz * (theta - theta_p) / theta + _radius * sin(theta_p) + h * sin(theta_p);
  }

  return disp;
}
