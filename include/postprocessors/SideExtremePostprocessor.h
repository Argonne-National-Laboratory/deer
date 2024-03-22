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
  virtual PostprocessorValue getValue() const;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

protected:
  unsigned int _qp;
  PostprocessorValue _curr_value;
  const VariableValue & _u;
  const ExtremeType _type;
};
