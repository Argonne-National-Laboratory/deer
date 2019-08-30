#ifndef THICKNESSGRADIENT_H
#define THICKNESSGRADIENT_H

#include "Function.h"

class ThicknessGradient;

template <>
InputParameters validParams<ThicknessGradient>();

class ThicknessGradient : public Function
{
 public:
  ThicknessGradient(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p);

 private:
  Real _getTemp(Real T2p, Real nx);

 protected:
  Real _delay, _T1, _T2, _tramp1, _thold1, _tramp2, _thold2, _x1, _x2;
  int _index;

};



#endif // THICKNESSGRADIENT_H
