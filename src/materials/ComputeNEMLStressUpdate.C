#include "ComputeNEMLStressUpdate.h"

registerMooseObject("DeerApp", ComputeNEMLStressUpdate);

template <>
InputParameters
validParams<ComputeNEMLStressUpdate>()
{
  InputParameters params = validParams<ComputeNEMLStressBase>();

  return params;
}


// Note: I am not really using use_displaced_mesh, the stress update
// model doesn't care where it's evaluated.
// use_displaced_mesh is only used to say "please use a large deformation
// stress update", consistent with the kernel and the strain material
ComputeNEMLStressUpdate::ComputeNEMLStressUpdate(
    const InputParameters & parameters)
  : ComputeNEMLStressBase(parameters),
    _ld(getParam<bool>("use_displaced_mesh"))
{

}

void
ComputeNEMLStressUpdate::stressUpdate(
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
  if (_ld) {
    ier = _model->update_ld_inc(e_np1, e_n, w_np1, w_n,
                                T_np1, T_n, t_np1, t_n, 
                                s_np1, s_n, h_np1, h_n, 
                                A_np1, B_np1, u_np1, u_n,
                                p_np1, p_n);
  }
  else {
    ier = _model->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n,
                      s_np1, s_n, h_np1, h_n, A_np1, u_np1, u_n,
                      p_np1, p_n);
    std::fill(B_np1, B_np1+18, 0.0);
  }

  if (ier != neml::SUCCESS)
    throw MooseException("NEML stress update failed!");
}

