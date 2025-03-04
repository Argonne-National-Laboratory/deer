#include "NEMLCrystalPlasticityFromAux.h"

registerMooseObject("DeerApp", NEMLCrystalPlasticityFromAux);

InputParameters
NEMLCrystalPlasticityFromAux::validParams()
{
  InputParameters params = CauchyStressFromNEML::validParams();

  params.addRequiredCoupledVar("initial_orientation", "Aux with initial orientations");

  return params;
}

NEMLCrystalPlasticityFromAux::NEMLCrystalPlasticityFromAux(const InputParameters & parameters)
  : CauchyStressFromNEML(parameters),
    _orientation(declareProperty<std::vector<Real>>("orientation")),
    _initial_orientation(coupledVectorValue("initial_orientation"))
{
  _cpmodel = static_cast<neml::SingleCrystalModel *>(_model.get());
}

void
NEMLCrystalPlasticityFromAux::initQpStatefulProperties()
{
  CauchyStressFromNEML::initQpStatefulProperties();

  auto q = neml::Orientation::createRodrigues(&(_initial_orientation[_qp](0)));
  _cpmodel->set_active_orientation(&_history[_qp].front(), q);
}

void
NEMLCrystalPlasticityFromAux::computeQpCauchyStress()
{
  CauchyStressFromNEML::computeQpCauchyStress();
  _formCPOutput();
}

void
NEMLCrystalPlasticityFromAux::_formCPOutput()
{
  _orientation[_qp].resize(4);
  neml::Orientation q = _cpmodel->get_active_orientation(&_history[_qp].front());
  std::copy(q.quat(), q.quat() + 4, _orientation[_qp].begin());
}
