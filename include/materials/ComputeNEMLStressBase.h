#pragma once

#include "ComputeNEMLStress.h"
#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

#include "neml_interface.h"

class ComputeNEMLStressBase : public DerivativeMaterialInterface<Material> {
public:
  static InputParameters validParams();
  ComputeNEMLStressBase(const InputParameters &parameters);
  virtual ~ComputeNEMLStressBase(){};

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  virtual void stressUpdate(const double *const e_np1, const double *const e_n,
                            const double *const w_np1, const double *const w_n,
                            double T_np1, double T_n, double t_np1, double t_n,
                            double *const s_np1, const double *const s_n,
                            double *const h_np1, const double *const h_n,
                            double *const A_np1, double *const B_np1,
                            double &u_np1, double u_n, double &p_np1,
                            double p_n) = 0;

private:
  void updateStrain();
  // Extract the polarization stress from the provided scalar values
  // Ordering (per MOOSE) is 11 22 33 23 13 12
  RankTwoTensor setupPolarizationStress();

protected:
  FileName _fname;
  std::string _mname;
  std::unique_ptr<neml::NEMLModel> _model;

  const MaterialProperty<RankTwoTensor> &_mechanical_strain_inc;
  const MaterialProperty<RankTwoTensor> &_vorticity_inc;

  const VariableValue &_temperature; // Will default to zero
  const VariableValue &_temperature_old;

  MaterialProperty<RankTwoTensor> &_mechanical_strain;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain_old;

  MaterialProperty<RankTwoTensor> &_linear_rot;
  const MaterialProperty<RankTwoTensor> &_linear_rot_old;

  MaterialProperty<RankTwoTensor> &_base_stress;
  const MaterialProperty<RankTwoTensor> &_base_stress_old;

  MaterialProperty<RankTwoTensor> &_stress;
  const MaterialProperty<RankTwoTensor> &_stress_old;

  MaterialProperty<RankFourTensor> &_material_jacobian;

  MaterialProperty<std::vector<Real>> &_hist;
  const MaterialProperty<std::vector<Real>> &_hist_old;

  MaterialProperty<Real> &_energy;
  const MaterialProperty<Real> &_energy_old;
  MaterialProperty<Real> &_dissipation;
  const MaterialProperty<Real> &_dissipation_old;

  MaterialProperty<RankTwoTensor> &_elastic_strain;
  MaterialProperty<RankTwoTensor> &_inelastic_strain;
 
  const bool _ld;

  unsigned int _num_polarization;
  std::vector<const VariableValue*> _polarization_stress;
};
