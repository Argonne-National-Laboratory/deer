#pragma once

#include "ComputeThermalExpansionEigenstrainBase.h"

#include "neml_interface.h"

/**
 *  ComputeThermalExpansionEigenstrainNEML computes the thermal expansion
 *  strain from the instantaneous CTE provided by a NEML model
 */
class ComputeThermalExpansionEigenstrainNEML : public ComputeThermalExpansionEigenstrainBase
{
public:
  static InputParameters validParams();
  ComputeThermalExpansionEigenstrainNEML(const InputParameters & parameters);
  virtual void initQpStatefulProperties() override;

protected:
  virtual void computeThermalStrain(Real & thermal_strain, Real & instantaneous_cte) override;

protected:
  FileName _fname;
  std::string _mname;
  std::unique_ptr<neml::NEMLModel> _model;

  MaterialProperty<Real> & _tstrain;
  const MaterialProperty<Real> & _tstrain_old;

  const VariableValue & _temperature_old;
};
