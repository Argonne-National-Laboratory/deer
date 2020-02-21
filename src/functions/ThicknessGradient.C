#include "ThicknessGradient.h"

registerMooseObject("DeerApp", ThicknessGradient);

InputParameters ThicknessGradient::validParams() {
  InputParameters params = Function::validParams();
  params.addParam<Real>("delay", 0.0, "Keep T1 until...");
  params.addParam<Real>("T1", 0.0, "Initial temperature");
  params.addParam<Real>("T2", 0.0, "Peak temperature");
  params.addParam<Real>("x1", 0.0, "Starting coordinate");
  params.addParam<Real>("x2", 0.0, "Ending coordinate");
  params.addParam<Real>("tramp1", 0.0, "Time to ramp up");
  params.addParam<Real>("thold1", 0.0, "Time to hold at high temperature");
  params.addParam<Real>("tramp2", 0.0, "Time to ramp down");
  params.addParam<Real>("thold2", 0.0, "Time to hold at low temperature");
  params.addParam<int>("index", 0, "Which direction to draw the gradient");

  return params;
}

ThicknessGradient::ThicknessGradient(const InputParameters &parameters)
    : Function(parameters), _delay(getParam<Real>("delay")),
      _T1(getParam<Real>("T1")), _T2(getParam<Real>("T2")),
      _tramp1(getParam<Real>("tramp1")), _thold1(getParam<Real>("thold1")),
      _tramp2(getParam<Real>("tramp2")), _thold2(getParam<Real>("thold2")),
      _x1(getParam<Real>("x1")), _x2(getParam<Real>("x2")),
      _index(getParam<int>("index")) {}

Real ThicknessGradient::value(Real t, const Point &p) {
  Real x = p(_index);
  Real period = _tramp1 + _thold1 + _tramp2 + _thold2;

  Real nt = fmod(t - _delay, period);
  Real nx = (x - _x1) / (_x2 - _x1);

  if (t > _delay) {
    Real T2p;
    if (nt < _tramp1) {
      T2p = _T1 + (_T2 - _T1) / _tramp1 * nt;
    } else if (nt < (_tramp1 + _thold1)) {
      T2p = _T2;
    } else if (nt < (_tramp1 + _thold1 + _tramp2)) {
      T2p = _T2 + (_T1 - _T2) / _tramp2 * (nt - _tramp1 - _thold1);
    } else {
      T2p = _T1;
    }

    return _getTemp(T2p, nx);
  } else {
    return _T1;
  }
}

Real ThicknessGradient::_getTemp(Real T2p, Real nx) {
  return (T2p - _T1) * nx + _T1;
}
