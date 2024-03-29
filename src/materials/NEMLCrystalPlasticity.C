#include "NEMLCrystalPlasticity.h"

registerMooseObject("DeerApp", NEMLCrystalPlasticity);

InputParameters
NEMLCrystalPlasticity::validParams()
{
  InputParameters params = CauchyStressFromNEML::validParams();

  params.addParam<UserObjectName>("euler_angle_reader", "Name of Euler angle reader");
  params.addParam<unsigned int>("grain_id", "ID of the grain for this material");
  params.addParam<std::string>("angle_type", "degrees", "Euler angle type: degrees or radians");
  params.addParam<std::string>(
      "angle_convention", "kocks", "Euler angle convention: kocks, bunge, or roe");
  return params;
}

NEMLCrystalPlasticity::NEMLCrystalPlasticity(const InputParameters & parameters)
  : CauchyStressFromNEML(parameters),
    _orientation(declareProperty<std::vector<Real>>("orientation")),
    _using_reader(parameters.isParamSetByUser("euler_angle_reader")),
    _euler_angles(_using_reader ? &getUserObjectByName<ElementPropertyReadFile>(
                                      getParam<UserObjectName>("euler_angle_reader"))
                                : nullptr),
    _angle_type(getParam<std::string>("angle_type")),
    _angle_convention(getParam<std::string>("angle_convention"))
{
  if (!_using_reader)
    mooseWarning("No Euler angle reader provided, using angles from NEML");

  _cpmodel = static_cast<neml::SingleCrystalModel *>(_model.get());
}

void
NEMLCrystalPlasticity::initQpStatefulProperties()
{
  CauchyStressFromNEML::initQpStatefulProperties();

  if (_using_reader)
  {
    auto q = neml::Orientation::createEulerAngles(_euler_angles->getData(_current_elem, 0),
                                                  _euler_angles->getData(_current_elem, 1),
                                                  _euler_angles->getData(_current_elem, 2),
                                                  _angle_type,
                                                  _angle_convention);
    _cpmodel->set_active_orientation(&_history[_qp].front(), q);
  }
}

void
NEMLCrystalPlasticity::computeQpCauchyStress()
{
  CauchyStressFromNEML::computeQpCauchyStress();
  _formCPOutput();
}

void
NEMLCrystalPlasticity::_formCPOutput()
{
  _orientation[_qp].resize(4);
  neml::Orientation q = _cpmodel->get_active_orientation(&_history[_qp].front());
  std::copy(q.quat(), q.quat() + 4, _orientation[_qp].begin());
}
