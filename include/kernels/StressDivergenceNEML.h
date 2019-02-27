#ifndef STRESSDIVERGENCENEML_H
#define STRESSDIVERGENCENEML_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class StressDivergenceNEML;

template <>
InputParameters validParams<StressDivergenceNEML>();

class StressDivergenceNEML: public DerivativeMaterialInterface<Kernel>
{
 public:
  StressDivergenceNEML(const InputParameters & parameters);
  virtual ~StressDivergenceNEML() {};
  
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
  Real matJacobianComponent(const RankFourTensor & A, const RankFourTensor & B,
                           const RankFourTensor & E, const RankFourTensor & W,
                           unsigned int i, unsigned int m,
                           const RealGradient & grad_psi,
                           const RealGradient & grad_phi);
  Real geomJacobianComponent(unsigned int i, unsigned int m,
                           const RealGradient & grad_psi,
                           const RealGradient & grad_phi,
                           const RankTwoTensor & stress);

 protected:
  bool _ld;

  unsigned int _component;
  unsigned int _ndisp;

  std::vector<unsigned int> _disp_nums;
  std::vector<MooseVariable*> _disp_vars;
  std::vector<const VariableGradient *> _grad_disp;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankFourTensor> & _material_strain_jacobian;
  const MaterialProperty<RankFourTensor> & _material_vorticity_jacobian;
  const MaterialProperty<RankFourTensor> & _strain_grad;
  const MaterialProperty<RankFourTensor> & _vorticity_grad;

  RankTwoTensor _qp_def_grad;
};

#endif // STRESSDIVERGENCENEML_H
