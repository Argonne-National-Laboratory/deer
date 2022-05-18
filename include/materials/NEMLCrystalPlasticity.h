#pragma once

#include "CauchyStressFromNEML.h"

#include "ElementPropertyReadFile.h"

#include "cp/singlecrystal.h"

class NEMLCrystalPlasticity : public CauchyStressFromNEML
{
public:
  static InputParameters validParams();
  NEMLCrystalPlasticity(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpCauchyStress();

private:
  void _formCPOutput();

protected:
  MaterialProperty<std::vector<Real>> & _orientation;
  const bool _using_reader;
  const ElementPropertyReadFile * _euler_angles;
  std::string _angle_type;
  std::string _angle_convention;
  neml::SingleCrystalModel * _cpmodel;
};
