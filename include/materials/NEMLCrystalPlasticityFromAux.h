#pragma once

#include "CauchyStressFromNEML.h"

#include "cp/singlecrystal.h"

class NEMLCrystalPlasticityFromAux : public CauchyStressFromNEML
{
public:
  static InputParameters validParams();
  NEMLCrystalPlasticityFromAux(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpCauchyStress();

private:
  void _formCPOutput();

protected:
  MaterialProperty<std::vector<Real>> & _orientation;
  const VectorVariableValue & _initial_orientation;
  neml::SingleCrystalModel * _cpmodel;
};
