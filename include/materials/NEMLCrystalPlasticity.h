#pragma once

#include "CauchyStressFromNEML.h"

#include "EulerAngleProvider.h"

#include "cp/singlecrystal.h"

class NEMLCrystalPlasticity: public CauchyStressFromNEML
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
   const EulerAngleProvider * _euler_angles;
   const bool _set_grain_id;
   unsigned int _grain_id;
   std::string _angle_type;
   std::string _angle_convention;
   neml::SingleCrystalModel * _cpmodel;
};
