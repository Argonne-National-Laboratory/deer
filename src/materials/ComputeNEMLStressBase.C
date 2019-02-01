#include "ComputeNEMLStressBase.h"

template <>
InputParameters
validParams<ComputeNEMLStressBase>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<FileName>("database", "Path to NEML XML database.");
  params.addRequiredParam<std::string>("model", "Model name in NEML database.");
  params.addCoupledVar("temperature", 0.0, "Coupled temperature");
  return params;
}

ComputeNEMLStressBase::ComputeNEMLStressBase(const InputParameters & parameters) 
  : DerivativeMaterialInterface<Material>(parameters),
    _fname(getParam<FileName>("database")),
    _mname(getParam<std::string>("model")),
    _mechanical_strain_inc(getMaterialPropertyByName<RankTwoTensor>("mechanical_strain_inc")),
    _vorticity_inc(getMaterialPropertyByName<RankTwoTensor>("vorticity_inc")),
    _temperature(coupledValue("temperature")),
    _temperature_old(coupledValueOld("temperature")),
    _mechanical_strain(declareProperty<RankTwoTensor>("mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>("mechanical_strain")),
    _stress(declareProperty<RankTwoTensor>("stress")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>("stress")),
    _material_strain_jacobian(declareProperty<RankFourTensor>("material_strain_jacobian")),
    _material_vorticity_jacobian(declareProperty<RankFourTensor>("material_vorticity_jacobian")),
    _hist(declareProperty<std::vector<Real>>("hist")),
    _hist_old(getMaterialPropertyOld<std::vector<Real>>("hist")),
    _energy(declareProperty<Real>("energy")),
    _energy_old(getMaterialPropertyOld<Real>("energy")),
    _dissipation(declareProperty<Real>("dissipation")),
    _dissipation_old(getMaterialPropertyOld<Real>("dissipation")),
    _shear_modulus(declareProperty<Real>("shear_modulus")),
    _bulk_modulus(declareProperty<Real>("bulk_modulus")),
    _elasticity_tensor(declareProperty<RankFourTensor>("elasticity_tensor")),
    _elastic_strain(declareProperty<RankTwoTensor>("elastic_strain")),
    _inelastic_strain(declareProperty<RankTwoTensor>("inelastic_strain"))
{
  // I strongly hesitate to put this here, may change later
  _model = neml::parse_xml_unique(_fname, _mname);
}

void ComputeNEMLStressBase::computeQpProperties()
{
  stressUpdate();
}

void ComputeNEMLStressBase::initQpStatefulProperties()
{
  // Basic variables maintained here
  _mechanical_strain[_qp].zero();
  _stress[_qp].zero();

  // Figure out initial history
  int ier;
  _hist[_qp].resize(_model->nhist());
  ier = _model->init_hist(&(_hist[_qp][0]));
  
  if (ier != neml::SUCCESS) {
    mooseError("Error initializing NEML history!");
  }

  // Various other junk
  _energy[_qp] = 0.0;
  _dissipation[_qp] = 0.0;
}
