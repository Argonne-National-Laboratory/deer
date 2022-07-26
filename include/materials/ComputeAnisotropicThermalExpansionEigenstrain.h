#pragma once

#include "Material.h"
#include "RankTwoTensorForward.h"
#include "Function.h"

class ComputeAnisotropicThermalExpansionEigenstrain : public Material
{
public:
  static InputParameters validParams();

  ComputeAnisotropicThermalExpansionEigenstrain(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeQpProperties() override;

  /// Base name prepended to material property name
  const std::string _base_name;

  /// Stores the current total eigenstrain
  MaterialProperty<RankTwoTensor> & _eigenstrain;

  const VariableValue & _T;
  const VariableValue & _T_old;

  const Function & _CTE0;
  const Function & _CTE1;
  const Function & _CTE2;

  const MaterialProperty<RealVectorValue> & _n0;
  const MaterialProperty<RealVectorValue> & _n1;
  const MaterialProperty<RealVectorValue> & _n2;

  MaterialProperty<Real> & _thermal_strain_0;
  MaterialProperty<Real> & _thermal_strain_1;
  MaterialProperty<Real> & _thermal_strain_2;

  const MaterialProperty<Real> & _thermal_strain_0_old;
  const MaterialProperty<Real> & _thermal_strain_1_old;
  const MaterialProperty<Real> & _thermal_strain_2_old;

  /// Restartable data to check for the zeroth and first time steps for thermal calculations
  bool & _step_zero;
};
