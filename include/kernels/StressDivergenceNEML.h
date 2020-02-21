#pragma once

#include "DerivativeMaterialInterface.h"
#include "Kernel.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

class StressDivergenceNEML : public DerivativeMaterialInterface<Kernel> {
public:
  static InputParameters validParams();
  StressDivergenceNEML(const InputParameters &parameters);
  virtual ~StressDivergenceNEML(){};

protected:
  virtual void initialSetup() override;

  // These must be overwritten for Bbar-type stabilizations
  virtual void precalculateResidual() override;
  virtual void precalculateJacobian() override;
  virtual void precalculateOffDiagJacobian(unsigned int jvar) override;

  // These are the standard implementation functions
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  Real matJacobianComponent(const RankFourTensor &C, unsigned int i,
                            unsigned int m, const RealGradient &grad_psi,
                            const RealGradient &grad_phi,
                            const RankTwoTensor &df);
  Real geomJacobianComponent(unsigned int i, unsigned int m,
                             const RealGradient &grad_psi,
                             const RealGradient &grad_phi,
                             const RankTwoTensor &stress);

protected:
  bool _ld;

  unsigned int _component;
  unsigned int _ndisp;

  std::vector<unsigned int> _disp_nums;
  std::vector<MooseVariable *> _disp_vars;
  std::vector<const VariableGradient *> _grad_disp;

  const MaterialProperty<RankTwoTensor> &_stress;
  const MaterialProperty<RankFourTensor> &_material_jacobian;

  const MaterialProperty<RankTwoTensor> &_df;
};

Real det(const RankTwoTensor &T);
