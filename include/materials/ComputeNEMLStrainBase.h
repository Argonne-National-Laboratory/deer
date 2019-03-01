#ifndef COMPUTENEMLSTRAINBASE_H
#define COMPUTENEMLSTRAINBASE_H

#include "Material.h"
#include "DerivativeMaterialInterface.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class ComputeNEMLStrainBase;

template <>
InputParameters validParams<ComputeNEMLStrainBase>();

class ComputeNEMLStrainBase: public DerivativeMaterialInterface<Material>
{
 public:
  ComputeNEMLStrainBase(const InputParameters & parameters);
  virtual ~ComputeNEMLStrainBase() {};

 protected:
  virtual void initialSetup() override;
  virtual void initQpStatefulProperties() override;
  virtual void computeProperties() override;
  virtual void computeQpStatefulProperties();

  virtual void precalculate();

 protected:
  unsigned int _ndisp;

  std::vector<const VariableValue *> _disp;
  std::vector<const VariableGradient *> _grad_disp;
  std::vector<const VariableValue *> _disp_old;
  std::vector<const VariableGradient *> _grad_disp_old;

  MaterialProperty<RankTwoTensor> & _strain_inc;
  MaterialProperty<RankTwoTensor> & _mechanical_strain_inc;
  MaterialProperty<RankTwoTensor> & _vorticity_inc;

  MaterialProperty<RankFourTensor> & _strain_grad;
  MaterialProperty<RankFourTensor> & _vorticity_grad;

  MaterialProperty<RankTwoTensor> & _def_grad_inv;
  const MaterialProperty<RankTwoTensor> & _def_grad_inv_old;

  bool _ld;
};

#endif // COMPUTENEMLSTRAINBASE_H
