#include "SideExtremeMaterialProperty.h"

registerMooseObject("DeerApp", SideExtremeMaterialProperty);

InputParameters
SideExtremeMaterialProperty::validParams()
{
  MooseEnum type_options("max=0 min=1", "max");

  InputParameters params = SidePostprocessor::validParams();

  params.addRequiredParam<MaterialPropertyName>("mat_prop",
                                                "Material property for which to find the extreme");
  params.addParam<MooseEnum>("value_type",
                             type_options,
                             "Type of extreme value to return. "
                             "Options are: 'max' and 'min'");

  params.addClassDescription("Finds the maximum value of variable over the side set");

  return params;
}

SideExtremeMaterialProperty::SideExtremeMaterialProperty(const InputParameters & parameters)
  : SidePostprocessor(parameters),
    _mat_prop(getMaterialProperty<Real>("mat_prop")),
    _qp(0),
    _curr_value(0.0),
    _type((ExtremeType)(int)parameters.get<MooseEnum>("value_type"))
{
}

void
SideExtremeMaterialProperty::initialize()
{
  switch (_type)
  {
    case MAX:
      _curr_value = -std::numeric_limits<Real>::infinity();
      break;
    case MIN:
      _curr_value = std::numeric_limits<Real>::infinity();
      break;
  }
}

void
SideExtremeMaterialProperty::execute()
{
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    switch (_type)
    {
      case MAX:
        _curr_value = std::max(_curr_value, _mat_prop[_qp]);
        break;
      case MIN:
        _curr_value = std::min(_curr_value, _mat_prop[_qp]);
        break;
    }
  }
}

Real
SideExtremeMaterialProperty::getValue()
{
  switch (_type)
  {
    case MAX:
      gatherMax(_curr_value);
      break;
    case MIN:
      gatherMin(_curr_value);
      break;
  }
  return _curr_value;
}

void
SideExtremeMaterialProperty::threadJoin(const UserObject & y)
{
  const SideExtremeMaterialProperty & pps = static_cast<const SideExtremeMaterialProperty &>(y);
  switch (_type)
  {
    case MAX:
      _curr_value = std::max(_curr_value, pps._curr_value);
      break;
    case MIN:
      _curr_value = std::min(_curr_value, pps._curr_value);
      break;
  }
}
