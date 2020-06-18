//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SalehaniIrani3DCTractionCZMIncremental.h"

registerMooseObject("DeerApp", SalehaniIrani3DCTractionCZMIncremental);

InputParameters SalehaniIrani3DCTractionCZMIncremental::validParams() {
  InputParameters params = CZMMaterialBasePK1::validParams();
  params.addClassDescription(
      "3D Coupled (3DC) cohesive law of Salehani and Irani with no damage");
  params.addRequiredParam<Real>(
      "normal_gap_at_maximum_normal_traction",
      "The value of normal gap at which maximum normal traction is achieved");
  params.addRequiredParam<Real>("tangential_gap_at_maximum_shear_traction",
                                "The value of tangential gap at which maximum "
                                "shear traction is achieved");
  params.addRequiredParam<Real>(
      "maximum_normal_traction",
      "The maximum normal traction the interface can sustain");
  params.addRequiredParam<Real>(
      "maximum_shear_traction",
      "The maximum shear traction the interface can sustain");
  return params;
}

SalehaniIrani3DCTractionCZMIncremental::SalehaniIrani3DCTractionCZMIncremental(
    const InputParameters &parameters)
    : CZMMaterialBasePK1(parameters),
      _delta_u0(
          {getParam<Real>("normal_gap_at_maximum_normal_traction"),
           std::sqrt(2) *
               getParam<Real>("tangential_gap_at_maximum_shear_traction"),
           std::sqrt(2) *
               getParam<Real>("tangential_gap_at_maximum_shear_traction")}),
      _max_allowable_traction({getParam<Real>("maximum_normal_traction"),
                               getParam<Real>("maximum_shear_traction"),
                               getParam<Real>("maximum_shear_traction")}) {}

void SalehaniIrani3DCTractionCZMIncremental::
    computeTractionIncrementAndDerivatives() {
  _traction_inc[_qp] = computeTraction();
  _dtraction_djump[_qp] = computeTractionDerivatives();
}

RealVectorValue SalehaniIrani3DCTractionCZMIncremental::computeTraction() {

  RealVectorValue traction_local_inc;

  Real A = 0;
  Real B = 0;

  for (unsigned int i = 0; i < 3; i++) {
    Real alpha = (i == 0) ? 1. : 2.;
    A -= std::pow(_displacement_jump[_qp](i) / _delta_u0[i], alpha);
    B -= alpha * _displacement_jump_inc[_qp](i) / _delta_u0[i] *
         std::pow(_displacement_jump[_qp](i) / _delta_u0[i], alpha - 1.);
    ;
  }
  A = std::exp(A);

  for (unsigned int i = 0; i < 3; i++) {
    Real lambda = (i == 0) ? std::exp(1.) : std::sqrt(2. * std::exp(1.));

    traction_local_inc(i) =
        _max_allowable_traction[i] * lambda / _delta_u0[i] * A *
        (_displacement_jump_inc[_qp](i) + _displacement_jump[_qp](i) * B);
  }

  return traction_local_inc;
}

RankTwoTensor
SalehaniIrani3DCTractionCZMIncremental::computeTractionDerivatives() {
  RankTwoTensor dtraction_djumpinc;

  Real A = 0;
  Real B = 0;
  for (unsigned int i = 0; i < 3; i++) {
    Real alpha = (i == 0) ? 1. : 2.;
    A -= std::pow(_displacement_jump[_qp](i) / _delta_u0[i], alpha);
    B -= alpha * _displacement_jump_inc[_qp](i) / _delta_u0[i] *
         std::pow(_displacement_jump[_qp](i) / _delta_u0[i], alpha - 1.);
  }
  A = std::exp(A);

  for (unsigned int i = 0; i < 3; i++) {
    Real lambda = (i == 0) ? std::exp(1.) : std::sqrt(2. * std::exp(1.));
    for (unsigned int j = 0; j < 3; j++) {
      Real alpha = (j == 0) ? 1. : 2.;
      Real dA_djump =
          A * (-alpha / _delta_u0[j] *
               std::pow(_displacement_jump[_qp](j) / _delta_u0[j], alpha - 1.));
      Real dB_djump =
          -alpha / _delta_u0[j] *
          std::pow(_displacement_jump[_qp](j) / _delta_u0[j], alpha - 1.);

      if (_displacement_jump_inc[_qp](j) != 0)
        dB_djump -=
            (alpha * alpha - alpha) / (_delta_u0[j] * _delta_u0[j]) *
            _displacement_jump_inc[_qp](j) *
            std::pow(_displacement_jump[_qp](j) / _delta_u0[j], alpha - 2.);

      Real c = (i == j) ? 1. + B : 0;
      dtraction_djumpinc(i, j) =
          dA_djump * (_displacement_jump_inc[_qp](i) +
                      _displacement_jump[_qp](i) * B) +
          A * (c + _displacement_jump[_qp](i) * dB_djump);

      dtraction_djumpinc(i, j) *=
          _max_allowable_traction[i] * lambda / _delta_u0[i];
    }
  }
  return dtraction_djumpinc;
}
