#pragma once

#include "Function.h"
#include "FunctionInterface.h"

class InductionFunction: public Function, public FunctionInterface {
 public:
  static InputParameters validParams();
  InductionFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const;

 protected:
  Point _center;
  Point _axis;
  Real _radius;
  Real _thickness;
  Real _delta;
};
