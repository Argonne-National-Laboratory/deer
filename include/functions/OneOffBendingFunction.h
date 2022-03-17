#pragma once

#include "Function.h"
#include "FunctionInterface.h"

class OneOffBendingFunction : public Function, public FunctionInterface {
 public:
  static InputParameters validParams();
  OneOffBendingFunction(const InputParameters & parameters);

  Real value(Real t, const Point & pt) const override;

 protected:
  RealVectorValue displacements(Real t, const Point & pt) const;

 protected:
  unsigned int _component;
  const Real _length;
  const Real _radius;
  const Function & _theta;
};
