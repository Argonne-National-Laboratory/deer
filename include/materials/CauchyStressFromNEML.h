#pragma once

#include "ComputeLagrangianStressCauchy.h"

#include "neml_interface.h"

class CauchyStressFromNEML : public ComputeLagrangianStressCauchy
{
public:
  static InputParameters validParams();
  CauchyStressFromNEML(const InputParameters & parameters);

protected:
  virtual void computeQpCauchyStress();
  virtual void initQpStatefulProperties();

protected:
  FileName _fname;
  std::string _mname;
  std::unique_ptr<neml::NEMLModel> _model;

  const VariableValue & _temperature;
  const VariableValue & _temperature_old;

  MaterialProperty<std::vector<Real>> & _history;
  const MaterialProperty<std::vector<Real>> & _history_old;

  MaterialProperty<Real> & _energy;
  const MaterialProperty<Real> & _energy_old;

  MaterialProperty<Real> & _dissipation;
  const MaterialProperty<Real> & _dissipation_old;

  MaterialProperty<RankTwoTensor> & _linear_rotation;
  const MaterialProperty<RankTwoTensor> & _linear_rotation_old;

  const MaterialProperty<RankTwoTensor> & _cauchy_stress_old;

  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  MaterialProperty<RankTwoTensor> & _inelastic_strain;
  MaterialProperty<RankTwoTensor> & _elastic_strain;

  MaterialProperty<Real> & _dissipation_rate;
};

/// Tensor -> Mandel
void tensor_neml(const RankTwoTensor & in, double * const out);

/// Mandel -> tensor
void neml_tensor(const double * const in, RankTwoTensor & out);

/// Tangent -> tensor
void neml_tangent(const double * const in, RankFourTensor & out);

/// Tensor -> skew vector
void tensor_skew(const RankTwoTensor & in, double * const out);

/// Skew vector -> tensor
void skew_tensor(const double * const in, RankTwoTensor & out);

/// Skew + symmetric parts to full tangent
void
recombine_tangent(const double * const Dpart, const double * const Wpart, RankFourTensor & out);
