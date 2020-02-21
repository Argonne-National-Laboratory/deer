#pragma once

#include "ComputeNEMLStressBase.h"

class ComputeNEMLStressUpdate : public ComputeNEMLStressBase {
public:
  static InputParameters validParams();
  ComputeNEMLStressUpdate(const InputParameters &parameters);
  virtual ~ComputeNEMLStressUpdate(){};

protected:
  virtual void stressUpdate(const double *const e_np1, const double *const e_n,
                            const double *const w_np1, const double *const w_n,
                            double T_np1, double T_n, double t_np1, double t_n,
                            double *const s_np1, const double *const s_n,
                            double *const h_np1, const double *const h_n,
                            double *const A_np1, double *const B_np1,
                            double &u_np1, double u_n, double &p_np1,
                            double p_n);
};
