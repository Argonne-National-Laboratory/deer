#include "SideExtremePostprocessor.h"

registerMooseObject("DeerApp", SideExtremePostprocessor);

InputParameters
SideExtremePostprocessor::validParams()
{
  MooseEnum type_options("max=0 min=1", "max");

  InputParameters params = SidePostprocessor::validParams();

  params.addRequiredCoupledVar("variable", "The name of the variable");
  params.addParam<MooseEnum>("value_type",
                             type_options,
                             "Type of extreme value to return. "
                             "Options are: 'max' and 'min'");

  params.addClassDescription("Finds the maximum value of variable over the side set");

  return params;
}

SideExtremePostprocessor::SideExtremePostprocessor(const InputParameters & parameters)
  : SidePostprocessor(parameters),
    _qp(0),
    _curr_value(0.0),
    _u(coupledValue("variable")),
    _type((ExtremeType)(int)parameters.get<MooseEnum>("value_type"))
{
}

void
SideExtremePostprocessor::initialize()
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
SideExtremePostprocessor::execute()
{
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    switch (_type)
    {
      case MAX:
        _curr_value = std::max(_curr_value, _u[_qp]);
        break;
      case MIN:
        _curr_value = std::min(_curr_value, _u[_qp]);
        break;
    }
  }
}

PostprocessorValue
SideExtremePostprocessor::getValue() const
{
  return _curr_value;
}

void
SideExtremePostprocessor::finalize()
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
}

void
SideExtremePostprocessor::threadJoin(const UserObject & y)
{
  const SideExtremePostprocessor & pps = static_cast<const SideExtremePostprocessor &>(y);
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
