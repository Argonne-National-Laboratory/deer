#include "CapGradient.h"

registerMooseObject("DeerApp", CapGradient);

InputParameters CapGradient::validParams() {
  InputParameters params = Function::validParams();
  params.addParam<Real>("delay", 0.0, "Keep T1 until...");
  params.addParam<Real>("T1", 0.0, "Initial temperature");
  params.addParam<Real>("T2", 0.0, "Peak temperature");
  params.addParam<Real>("r1", 0.0, "Starting coordinate");
  params.addParam<Real>("r2", 0.0, "Ending coordinate");
  params.addParam<Real>("trans", 0.0, "Transition height to cap");
  params.addParam<Real>("radius", 0.0, "Spherical cap radius");
  params.addParam<Real>("tramp1", 0.0, "Time to ramp up");
  params.addParam<Real>("thold1", 0.0, "Time to hold at high temperature");
  params.addParam<Real>("tramp2", 0.0, "Time to ramp down");
  params.addParam<Real>("thold2", 0.0, "Time to hold at low temperature");
  params.addParam<int>("index", 2, "Height direction");

  return params;
}

CapGradient::CapGradient(const InputParameters &parameters)
    : Function(parameters), _delay(getParam<Real>("delay")),
      _T1(getParam<Real>("T1")), _T2(getParam<Real>("T2")),
      _tramp1(getParam<Real>("tramp1")), _thold1(getParam<Real>("thold1")),
      _tramp2(getParam<Real>("tramp2")), _thold2(getParam<Real>("thold2")),
      _r1(getParam<Real>("r1")), _r2(getParam<Real>("r2")),
      _trans(getParam<Real>("trans")), _radius(getParam<Real>("radius")),
      _index(getParam<int>("index")) {}

Real CapGradient::value(Real t, const Point &p) {
  // This trusts that z = 0.0 for 2d meshes...
  double z;
  int oc[2];
  if (_index == 0) {
    z = p(0);
    oc[0] = 1;
    oc[1] = 2;
  } else if (_index == 1) {
    z = p(1);
    oc[0] = 0;
    oc[1] = 2;
  } else if (_index == 2) {
    z = p(2);
    oc[0] = 0;
    oc[1] = 1;
  } else {
    mooseError("Coordinate index provided is > 2!");
  }

  double r;
  if (z < _trans) {
    r = 0.0;
    for (int i = 0; i < 2; i++)
      r += p(oc[i]) * p(oc[i]);
    r = sqrt(r);
  } else {
    r = 0.0;
    for (int i = 0; i < 2; i++)
      r += p(oc[i]) * p(oc[i]);
    r += pow(p(_index) - _trans, 2.0);
    r = sqrt(r);
  }

  Real period = _tramp1 + _thold1 + _tramp2 + _thold2;

  Real nt = fmod(t - _delay, period);
  Real nx = (r - _r1) / (_r2 - _r1);

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

Real CapGradient::_getTemp(Real T2p, Real nx) { return (T2p - _T1) * nx + _T1; }
