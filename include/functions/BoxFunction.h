#pragma once

#include "Function.h"

class BoxFunction : public Function
{
public:
  static InputParameters validParams();
  BoxFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const;

protected:
  Point _lb, _ub;
  Real _value;
};
