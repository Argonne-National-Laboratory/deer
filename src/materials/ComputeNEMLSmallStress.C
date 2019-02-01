#include "ComputeNEMLSmallStress.h"

registerMooseObject("DeerApp", ComputeNEMLSmallStress);

template <>
InputParameters
validParams<ComputeNEMLSmallStress>()
{
  InputParameters params = validParams<ComputeNEMLStressBase>();
  return params;
}

ComputeNEMLSmallStress::ComputeNEMLSmallStress(
    const InputParameters & parameters)
  : ComputeNEMLStressBase(parameters)
{

}

void
ComputeNEMLSmallStress::stressUpdate()
{
  // Get the strains added together
  updateStrain();

  // First do some declaration and translation
  double s_np1[6];
  double s_n[6];
  tensor_neml(_stress_old[_qp], s_n);

  double e_np1[6];
  tensor_neml(_mechanical_strain[_qp], e_np1);
  double e_n[6];
  tensor_neml(_mechanical_strain_old[_qp], e_n);

  double t_np1 = _t;
  double t_n = _t - _dt;

  double T_np1 = _temperature[_qp];
  double T_n = _temperature_old[_qp];
  
  double * h_np1 = &(_hist[_qp][0]);
  const double * const h_n = &(_hist_old[_qp][0]);

  double A_np1[36];
  
  double u_np1;
  double u_n = _energy_old[_qp];

  double p_np1;
  double p_n = _dissipation_old[_qp];

  int ier;

  // Actually call the update
  ier = _model->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n,
                    s_np1, s_n, h_np1, h_n, A_np1, u_np1, u_n,
                    p_np1, p_n);

  if (ier != neml::SUCCESS)
    throw MooseException("NEML stress update failed!");

  // Do more translation, now back to tensors
  neml_tensor(s_np1, _stress[_qp]);
  neml_tangent(A_np1, _material_strain_jacobian[_qp]);

  // The vorticity doesn't get used
  _material_vorticity_jacobian[_qp].zero();

  // Get the elastic strain
  double estrain[6];
  ier = _model->elastic_strains(s_np1, T_np1, h_np1, estrain);

  if (ier != neml::SUCCESS)
    mooseError("Error in NEML call for elastic strain!");

  // Translate
  neml_tensor(estrain, _elastic_strain[_qp]);

  // For EPP purposes calculate the inelastic strain
  double pstrain[6];
  for (int i=0; i<6; i++) {
    pstrain[i] = e_np1[i] - estrain[i];
  }
  neml_tensor(pstrain, _inelastic_strain[_qp]);

  // Store dissipation
  _energy[_qp] = u_np1;
  _dissipation[_qp] = p_np1;

  // Store elastic properties at current time
  double mu = _model->shear(T_np1);
  double K = _model->bulk(T_np1);
  double l = K - 2.0 * mu / 3.0;
  
  _shear_modulus[_qp] = mu;
  _bulk_modulus[_qp] = K;

  std::vector<Real> props({l, mu});
  _elasticity_tensor[_qp].fillFromInputVector(props, RankFourTensor::FillMethod::symmetric_isotropic);
}

void
ComputeNEMLSmallStress::updateStrain()
{
  _mechanical_strain[_qp] = _mechanical_strain_old[_qp] +
      _mechanical_strain_inc[_qp];
}
