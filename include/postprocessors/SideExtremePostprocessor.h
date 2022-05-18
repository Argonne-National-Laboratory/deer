#pragma once

#include "SidePostprocessor.h"

class SideExtremePostprocessor : public SidePostprocessor
{
public:
  static InputParameters validParams();

  enum ExtremeType
  {
    MAX,
    MIN
  };

  SideExtremePostprocessor(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  unsigned int _qp;
  Real _curr_value;
  const VariableValue & _u;
  const ExtremeType _type;
};
