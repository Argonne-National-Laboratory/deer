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
    _linear_rot(declareProperty<RankTwoTensor>("linear_rot")),
    _linear_rot_old(getMaterialPropertyOld<RankTwoTensor>("linear_rot")),
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

void ComputeNEMLStressBase::computeProperties()
{
  if (_bnd) return;
  DerivativeMaterialInterface<Material>::computeProperties();
}

void ComputeNEMLStressBase::computeQpProperties()
{
  // Get the strains added together
  updateStrain();

  // First do some declaration and translation
  double s_np1[6];
  double s_n[6];
  tensor_neml(_stress_old[_qp], s_n);

  double e_np1[6];
  tensor_neml(_mechanical_strain[_qp], e_np1);
  double e_n[6];
  tensor_neml(_mechanical_strain_old[_qp], e_n);

  // vorticity
  double w_np1[3];
  tensor_skew(_linear_rot[_qp], w_np1);
  double w_n[3];
  tensor_skew(_linear_rot[_qp], w_n);

  double t_np1 = _t;
  double t_n = _t - _dt;

  double T_np1 = _temperature[_qp];
  double T_n = _temperature_old[_qp];
  
  double * h_np1;
  const double * h_n;
  
  // MOOSE vector debug doesn't like this
  if (_model->nstore() > 0) {
    h_np1 = &(_hist[_qp][0]);
    h_n = &(_hist_old[_qp][0]);
  }
  else {
    h_np1 = nullptr;
    h_n = nullptr;
  }

  double A_np1[36];
  double B_np1[18];
  
  double u_np1;
  double u_n = _energy_old[_qp];

  double p_np1;
  double p_n = _dissipation_old[_qp];

  stressUpdate(e_np1, e_n, w_np1, w_n, T_np1, T_n, t_np1, t_n,
               s_np1, s_n, h_np1, h_n, A_np1, B_np1, u_np1, u_n,
               p_np1, p_n);

  // Do more translation, now back to tensors
  neml_tensor(s_np1, _stress[_qp]);
  neml_tangent(A_np1, _material_strain_jacobian[_qp]);
  neml_skew_tangent(B_np1, _material_vorticity_jacobian[_qp]);

  // Get the elastic strain
  double estrain[6];
  int ier = _model->elastic_strains(s_np1, T_np1, h_np1, estrain);

  if (ier != neml::SUCCESS)
    mooseError("Error in NEML call for elastic strain!");

  // Translate
  neml_tensor(estrain, _elastic_strain[_qp]);

  // For EPP purposes calculate the inelastic strain
  double pstrain[6];
  for (int i=0; i<6; i++) {
    pstrain[i] = e_np1[i] - estrain[i];
  }
  neml_tensor(pstrain, _inelastic_strain[_qp]);

  // Store dissipation
  _energy[_qp] = u_np1;
  _dissipation[_qp] = p_np1;

  // Store elastic properties at current time
  double mu = _model->shear(T_np1);
  double K = _model->bulk(T_np1);
  double l = K - 2.0 * mu / 3.0;
  
  _shear_modulus[_qp] = mu;
  _bulk_modulus[_qp] = K;

  std::vector<Real> props({l, mu});
  _elasticity_tensor[_qp].fillFromInputVector(props, RankFourTensor::FillMethod::symmetric_isotropic);
}

void ComputeNEMLStressBase::initQpStatefulProperties()
{
  // Basic variables maintained here
  _mechanical_strain[_qp].zero();
  _linear_rot[_qp].zero();
  _stress[_qp].zero();

  // Figure out initial history
  int ier;
  _hist[_qp].resize(_model->nstore());
  if (_model->nstore() > 0) {
    ier = _model->init_store(&_hist[_qp].front());
  }
  else {
    ier = 0;
  }

  if (ier != neml::SUCCESS) {
    mooseError("Error initializing NEML history!");
  }

  // Various other junk
  _energy[_qp] = 0.0;
  _dissipation[_qp] = 0.0;
}

void
ComputeNEMLStressBase::updateStrain()
{
  _mechanical_strain[_qp] = _mechanical_strain_old[_qp] +
      _mechanical_strain_inc[_qp];
  _linear_rot[_qp] = _linear_rot_old[_qp] + 
      _vorticity_inc[_qp];
}
