#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

#include "neml_interface.h"

class ComputeNEMLStress : public DerivativeMaterialInterface<Material> {
public:
  static InputParameters validParams();
  ComputeNEMLStress(const InputParameters &parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

protected:
  MaterialProperty<RankTwoTensor> &_stress;
  MaterialProperty<RankTwoTensor> &_elastic_strain;

  const MaterialProperty<RankTwoTensor> &_mechanical_strain;
  const MaterialProperty<RankTwoTensor> &_extra_stress;

  MaterialProperty<RankFourTensor> &_Jacobian_mult;

  FileName _fname;
  std::string _mname;
  std::unique_ptr<neml::NEMLModel> _model;
  MaterialProperty<std::vector<Real>> &_hist;
  const MaterialProperty<std::vector<Real>> &_hist_old;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain_old;
  const MaterialProperty<RankTwoTensor> &_stress_old;
  MaterialProperty<Real> &_energy;
  const MaterialProperty<Real> &_energy_old;
  MaterialProperty<Real> &_dissipation;
  const MaterialProperty<Real> &_dissipation_old;
  const VariableValue &_temperature; // Will default to zero
  const VariableValue &_temperature_old;
  MaterialProperty<RankTwoTensor> &_inelastic_strain;
};

/// Tensor -> my notation
void tensor_neml(const RankTwoTensor &in, double *const out);

/// Vector -> tensor
void neml_tensor(const double *const in, RankTwoTensor &out);

/// Tangent -> tensor
void neml_tangent(const double *const in, RankFourTensor &out);

/// Tensor -> skew vector
void tensor_skew(const RankTwoTensor &in, double *const out);

/// Skew vector -> tensor
void skew_tensor(const double *const in, RankTwoTensor &out);

/// Skew + symmetric parts to full tangent
void recombine_tangent(const double *const Dpart, const double *const Wpart,
                       RankFourTensor &out);
