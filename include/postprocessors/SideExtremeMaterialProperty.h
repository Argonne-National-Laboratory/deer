#pragma once

#include "SidePostprocessor.h"

class SideExtremeMaterialProperty : public SidePostprocessor
{
public:
  static InputParameters validParams();

  enum ExtremeType
  {
    MAX,
    MIN
  };

  SideExtremeMaterialProperty(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  const MaterialProperty<Real> & _mat_prop;
  unsigned int _qp;
  Real _curr_value;
  const ExtremeType _type;
};
