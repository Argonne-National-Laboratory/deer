#include "ComputeNEMLSmallStress.h"

registerMooseObject("DeerApp", ComputeNEMLSmallStress);

template <>
InputParameters
validParams<ComputeNEMLSmallStress>()
{
  InputParameters params = validParams<ComputeNEMLStressBase>();

  params.set<bool>("use_displaced_mesh") = false;  

  return params;
}

ComputeNEMLSmallStress::ComputeNEMLSmallStress(
    const InputParameters & parameters)
  : ComputeNEMLStressBase(parameters)
{

}

void
ComputeNEMLSmallStress::stressUpdate(
      const double * const e_np1, const double * const e_n,
      const double * const w_np1, const double * const w_n,
      double T_np1, double T_n, double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1, double * const B_np1,
      double & u_np1, double u_n, double & p_np1, double p_n)
{
  int ier;

  // Actually call the update
  ier = _model->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n,
                    s_np1, s_n, h_np1, h_n, A_np1, u_np1, u_n,
                    p_np1, p_n);

  std::fill(B_np1, B_np1+18, 0.0);

  if (ier != neml::SUCCESS)
    throw MooseException("NEML stress update failed!");
}

