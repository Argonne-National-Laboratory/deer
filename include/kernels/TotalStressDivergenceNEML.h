#pragma once

#include "DerivativeMaterialInterface.h"
#include "Kernel.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

#include "HomogenizationConstraintIntegral.h" // Index constants
#include "MooseVariableScalar.h"

class TotalStressDivergenceNEML : public DerivativeMaterialInterface<Kernel> {
public:
  static InputParameters validParams();
  TotalStressDivergenceNEML(const InputParameters &parameters);
  virtual ~TotalStressDivergenceNEML(){};

protected:
  virtual void initialSetup() override;

  // These are the standard implementation functions
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;

private:
  Real largeDeformationResidual(const RealGradient & grad_phi);
  Real smallDeformationResidual(const RealGradient & grad_phi);
  
  Real smallDeformationMatJac(unsigned int i, unsigned int k,
                              const RealGradient & grad_phi,
                              const RealGradient & grad_psi);
  Real largeDeformationMatJac(unsigned int i, unsigned int k,
                              const RealGradient & grad_phi,
                              const RealGradient & grad_psi);
  Real largeDeformationGeoJac(unsigned int i, unsigned int k,
                              const RealGradient & grad_phi,
                              const RealGradient & grad_psi);
  Real computeBaseJacobian();
  Real computeConstraintJacobian();

  Real sdBaseJacobian();
  Real ldBaseJacobian();

  Real sdConstraintJacobianStrain();
  Real sdConstraintJacobianStress();
  Real ldConstraintJacobianStrain();
  Real ldConstraintJacobianStress();

protected:
  bool _ld;

  unsigned int _component;
  unsigned int _ndisp;

  std::vector<unsigned int> _disp_nums;
  std::vector<MooseVariable *> _disp_vars;
  std::vector<const VariableGradient *> _grad_disp;

  const MaterialProperty<RankTwoTensor> &_stress;
  const MaterialProperty<RankFourTensor> &_material_jacobian;

  const MaterialProperty<RankTwoTensor> &_inv_def_grad;
  const MaterialProperty<Real> &_detJ;
  const MaterialProperty<RankTwoTensor> &_df;

  unsigned int _macro_gradient_num;
  const MooseVariableScalar * _macro_gradient;
  const HomogenizationConstants::index_list _indices;
  
  std::vector<HomogenizationConstants::ConstraintType> _ctypes;

  unsigned int _h;
};
