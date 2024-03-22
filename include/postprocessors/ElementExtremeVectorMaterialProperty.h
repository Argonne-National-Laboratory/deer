#pragma once

#include "ElementPostprocessor.h"

class ElementExtremeVectorMaterialProperty : public ElementPostprocessor
{
public:
  static InputParameters validParams();

  enum ExtremeType
  {
    MAX,
    MIN
  };

  ElementExtremeVectorMaterialProperty(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual PostprocessorValue getValue() const;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

protected:
  virtual void computeQpValue();

  const MaterialProperty<std::vector<Real>> & _mat_prop;
  unsigned int _index;

  ExtremeType _type;
  PostprocessorValue _value;
  unsigned int _qp;
};
