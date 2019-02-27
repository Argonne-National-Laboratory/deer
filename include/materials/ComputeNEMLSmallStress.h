#ifndef COMPUTENEMLSMALLSTRESS_H
#define COMPUTENEMLSMALLSTRESS_H

#include "ComputeNEMLStressBase.h"

class ComputeNEMLSmallStress;

template <>
InputParameters validParams<ComputeNEMLSmallStress>();

class ComputeNEMLSmallStress: public ComputeNEMLStressBase
{
 public:
  ComputeNEMLSmallStress(const InputParameters & parameters);
  virtual ~ComputeNEMLSmallStress() {};

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
