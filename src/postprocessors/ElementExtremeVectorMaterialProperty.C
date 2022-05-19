#include "ElementExtremeVectorMaterialProperty.h"

registerMooseObject("DeerApp", ElementExtremeVectorMaterialProperty);

InputParameters
ElementExtremeVectorMaterialProperty::validParams()
{
  InputParameters params = ElementPostprocessor::validParams();

  params.addRequiredParam<MaterialPropertyName>("mat_prop", "Material property");
  MooseEnum type_options("max=0 min=1");
  params.addRequiredParam<MooseEnum>("value_type", type_options, "max or min");

  params.addRequiredParam<unsigned int>("index", "The index to consider for this kernel");

  params.addClassDescription("Max or min of a vector material property");

  return params;
}

ElementExtremeVectorMaterialProperty::ElementExtremeVectorMaterialProperty(
    const InputParameters & parameters)
  : ElementPostprocessor(parameters),
    _mat_prop(getMaterialProperty<std::vector<Real>>("mat_prop")),
    _index(parameters.get<unsigned int>("index")),
    _type((ExtremeType)(int)parameters.get<MooseEnum>("value_type")),
    _value(_type == 0 ? -std::numeric_limits<Real>::max() : std::numeric_limits<Real>::max()),
    _qp(0)
{
}

void
ElementExtremeVectorMaterialProperty::initialize()
{
  switch (_type)
  {
    case MAX:
      _value = -std::numeric_limits<Real>::max(); // start w/ the min
      break;

    case MIN:
      _value = std::numeric_limits<Real>::max(); // start w/ the max
      break;
  }
}

void
ElementExtremeVectorMaterialProperty::execute()
{
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    computeQpValue();
}

void
ElementExtremeVectorMaterialProperty::computeQpValue()
{
  switch (_type)
  {
    case MAX:
      _value = std::max(_value, _mat_prop[_qp][_index]);
      break;

    case MIN:
      _value = std::min(_value, _mat_prop[_qp][_index]);
      break;
  }
}

Real
ElementExtremeVectorMaterialProperty::getValue()
{
  switch (_type)
  {
    case MAX:
      gatherMax(_value);
      break;
    case MIN:
      gatherMin(_value);
      break;
  }

  return _value;
}

void
ElementExtremeVectorMaterialProperty::threadJoin(const UserObject & y)
{
  const ElementExtremeVectorMaterialProperty & pps =
      static_cast<const ElementExtremeVectorMaterialProperty &>(y);

  switch (_type)
  {
    case MAX:
      _value = std::max(_value, pps._value);
      break;
    case MIN:
      _value = std::min(_value, pps._value);
      break;
  }
}
