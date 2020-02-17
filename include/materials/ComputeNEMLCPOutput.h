#pragma once

#include "ComputeNEMLStressUpdate.h"
#include "EulerAngleProvider.h"

class ComputeNEMLCPOutput : public ComputeNEMLStressUpdate {
public:
  static InputParameters validParams();
  ComputeNEMLCPOutput(const InputParameters &parameters);
  virtual ~ComputeNEMLCPOutput(){};

protected:
  neml::SingleCrystalModel *_cpmodel = nullptr;

  MaterialProperty<std::vector<Real>> &_orientation_q;
  /// object providing the Euler angles
  const EulerAngleProvider *_euler;
  /// grain id
  unsigned int _grain;
  unsigned int _given = 1;

  virtual void initQpStatefulProperties() override;

  virtual void stressUpdate(const double *const e_np1, const double *const e_n,
                            const double *const w_np1, const double *const w_n,
                            double T_np1, double T_n, double t_np1, double t_n,
                            double *const s_np1, const double *const s_n,
                            double *const h_np1, const double *const h_n,
                            double *const A_np1, double *const B_np1,
                            double &u_np1, double u_n, double &p_np1,
                            double p_n);

  void getCPOutput(double *const h_np1);
};
