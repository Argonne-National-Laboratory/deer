#include "ComputeRadiationSwellingEigenstrain.h"

#include "Function.h"
#include "RankTwoTensor.h"

template <>
InputParameters
validParams<ComputeRadiationSwellingEigenstrain>()
{
  InputParameters params = validParams<ComputeEigenstrainBase>();
  params.addRequiredParam<FunctionName>("swelling", "Swelling as a function of dose");
  params.addRequiredParam<FunctionName>("dose_rate", "Dose rate as a function of time");
  return params;
}

ComputeRadiationSwellingEigenstrain::ComputeRadiationSwellingEigenstrain(
    const InputParameters & parameters) :
    ComputeEigenstrainBase(parameters),
    _swelling(getFunction("swelling")),
    _dose_rate(getFunction("dose_rate")),
    _dose(declareProperty<Real>(_base_name + "dose")),
    _dose_old(declarePropertyOld<Real>(_base_name + "dose"))
{

}

void ComputeRadiationSwellingEigenstrain::initQpStatefulProperties()
{
  ComputeEigenstrainBase::initQpStatefulProperties();
  _dose[_qp] = 0.0;
}

void ComputeRadiationSwellingEigenstrain::computeQpEigenstrain()
{
  // Update dose
  _dose[_qp] = _dose_old[_qp] + _dose_rate.value(_t, _q_point[_qp]) * _dt;

  // Get swell
  double swell = _swelling.value(_dose[_qp], _q_point[_qp]);

  // Stick in the eigenstrain
  for (unsigned i = 0; i < LIBMESH_DIM; ++i) {
    for (unsigned j = 0; j < LIBMESH_DIM; ++j) {
      if (i == j) {
        _eigenstrain[_qp](i, j) = swell / 3.0;
      }
      else {
        _eigenstrain[_qp](i, j) = 0.0;
      }
    }
  }
}
