#ifndef COMPUTENEMLLARGESTRESS_H
#define COMPUTENEMLLARGESTRESS_H

#include "ComputeNEMLStressBase.h"

class ComputeNEMLLargeStress;

template <>
InputParameters validParams<ComputeNEMLLargeStress>();

class ComputeNEMLLargeStress: public ComputeNEMLStressBase
{
 public:
  ComputeNEMLLargeStress(const InputParameters & parameters);
  virtual ~ComputeNEMLLargeStress() {};

 protected:
  virtual void stressUpdate(
      const double * const e_np1, const double * const e_n,
      const double * const w_np1, const double * const w_n,
      double T_np1, double T_n, double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1, double * const B_np1,
      double & u_np1, double u_n, double & p_np1, double p_n);

};

#endif
