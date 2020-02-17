#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"

#include "RankFourTensor.h"
#include "RankTwoTensor.h"

class ComputeNEMLStrainBase : public DerivativeMaterialInterface<Material> {
public:
  static InputParameters validParams();
  ComputeNEMLStrainBase(const InputParameters &parameters);
  virtual ~ComputeNEMLStrainBase(){};

protected:
  virtual void initialSetup() override;
  virtual void initQpStatefulProperties() override;
  virtual void computeProperties() override;
  virtual void computeQpProperties() override;

  virtual void precalculate();

  RankTwoTensor eigenstrainIncrement();

protected:
  unsigned int _ndisp;

  std::vector<const VariableValue *> _disp;
  std::vector<const VariableGradient *> _grad_disp;
  std::vector<const VariableValue *> _disp_old;
  std::vector<const VariableGradient *> _grad_disp_old;

  MaterialProperty<RankTwoTensor> &_strain_inc;
  MaterialProperty<RankTwoTensor> &_mechanical_strain_inc;
  MaterialProperty<RankTwoTensor> &_vorticity_inc;

  MaterialProperty<RankTwoTensor> &_def_grad;
  const MaterialProperty<RankTwoTensor> &_def_grad_old;

  MaterialProperty<RankTwoTensor> &_df;

  std::vector<MaterialPropertyName> _eigenstrain_names;
  std::vector<const MaterialProperty<RankTwoTensor> *> _eigenstrains;
  std::vector<const MaterialProperty<RankTwoTensor> *> _eigenstrains_old;

  bool _ld;
};
