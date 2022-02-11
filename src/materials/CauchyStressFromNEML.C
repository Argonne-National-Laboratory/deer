#include "CauchyStressFromNEML.h"

registerMooseObject("DeerApp", CauchyStressFromNEML);

InputParameters
CauchyStressFromNEML::validParams() {
  InputParameters params = ComputeLagrangianStressCauchy::validParams();

  params.addRequiredParam<FileName>("database", "Path to NEML XML database.");
  params.addRequiredParam<std::string>("model", "Model name in NEML database.");
  params.addCoupledVar("temperature", 0.0, "Coupled temperature");

  return params;
}

CauchyStressFromNEML::CauchyStressFromNEML(const InputParameters & parameters)
  : ComputeLagrangianStressCauchy(parameters),
    _fname(getParam<FileName>("database")),
    _mname(getParam<std::string>("model")),
    _temperature(coupledValue("temperature")),
    _temperature_old(coupledValueOld("temperature")),
    _history(declareProperty<std::vector<Real>>(_base_name + "history")),
    _history_old(getMaterialPropertyOld<std::vector<Real>>(_base_name + "history")),
    _energy(declareProperty<Real>(_base_name + "energy")),
    _energy_old(getMaterialPropertyOld<Real>(_base_name + "energy")),
    _dissipation(declareProperty<Real>(_base_name + "dissipation")),
    _dissipation_old(declareProperty<Real>(_base_name + "dissipation_old")),
    _linear_rotation(declareProperty<RankTwoTensor>(_base_name + "linear_rotation")),
    _linear_rotation_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "linear_rotation")),
    _cauchy_stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name +
                                                             "cauchy_stress")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>(_base_name +
                                                                "mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + 
                                                                "mechanical_strain")),
    _inelastic_strain(declareProperty<RankTwoTensor>(_base_name + "inelastic_strain")),
    _elastic_strain(declareProperty<RankTwoTensor>(_base_name +
                                                   "elastic_strain"))
{
  // Should raise an exception if it does not work
  _model = neml::parse_xml_unique(_fname, _mname);
}

void
CauchyStressFromNEML::computeQpCauchyStress()
{
  // Setup all the Mandel notation things we need
  double s_np1[6];
  double s_n[6];
  tensor_neml(_cauchy_stress_old[_qp], s_n);

  // Strain
  double e_np1[6];
  tensor_neml(_mechanical_strain[_qp], e_np1);
  double e_n[6];
  tensor_neml(_mechanical_strain_old[_qp], e_n);
  
  // Vorticity
  RankTwoTensor L;
  if (_large_kinematics) {
    L = RankTwoTensor::Identity() - _inv_df[_qp];
  }
  else {
    L.zero();
  }
  _linear_rotation[_qp] = _linear_rotation_old[_qp] + (L - L.transpose()) / 2.0;

  double w_np1[3];
  tensor_skew(_linear_rotation[_qp], w_np1);
  double w_n[3];
  tensor_skew(_linear_rotation_old[_qp], w_n);
  
  // Time
  double t_np1 = _t;
  double t_n = _t - _dt;
  
  // Temperature
  double T_np1 = _temperature[_qp];
  double T_n = _temperature_old[_qp];
  
  // Internal state
  double * h_np1;
  const double * h_n;

  // Just to keep MOOSE debug happy
  if (_model->nstore() > 0) {
    h_np1 = &(_history[_qp][0]);
    h_n = &(_history_old[_qp][0]);
  }
  else {
    h_np1 = nullptr;
    h_n = nullptr;
  }
  
  // Energy
  double u_np1;
  double u_n = _energy_old[_qp];
  
  // Dissipation
  double p_np1;
  double p_n = _dissipation_old[_qp];

  // Tangent
  double A_np1[36];
  double B_np1[18];

  // Call NEML!
  int ier;
  if (_large_kinematics) {
    ier = _model->update_ld_inc(e_np1, e_n, w_np1, w_n, T_np1, T_n, t_np1, t_n,
                                s_np1, s_n, h_np1, h_n, A_np1, B_np1, u_np1, 
                                u_n, p_np1, p_n);
  }
  else {
    ier = _model->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n, 
                            h_np1, h_n, A_np1, u_np1, u_n, p_np1, p_n);
    std::fill(B_np1, B_np1 + 18, 0.0);
  }
  if (ier != neml::SUCCESS)
    throw MooseException("NEML stress update failed!");

  double estrain[6];
  ier = _model->elastic_strains(s_np1, T_np1, h_np1, estrain);
  if (ier != neml::SUCCESS)
    throw MooseException("Error in NEML call for elastic strains!");

  // Translate back from Mandel notation
  neml_tensor(s_np1, _cauchy_stress[_qp]);
  recombine_tangent(A_np1, B_np1, _cauchy_jacobian[_qp]);
  _energy[_qp] = u_np1;
  _dissipation[_qp] = p_np1;
  
  neml_tensor(estrain, _elastic_strain[_qp]);
  _inelastic_strain[_qp] = _mechanical_strain[_qp] - _elastic_strain[_qp];
}

void
CauchyStressFromNEML::initQpStatefulProperties()
{
  ComputeLagrangianStressCauchy::initQpStatefulProperties();
  
  int ier;
  _history[_qp].resize(_model->nstore());
  // This is only needed because MOOSE whines about zero sized vectors
  // that are not initialized
  if (_history[_qp].size() > 0)
    ier = _model->init_store(&_history[_qp].front());
  if (ier != neml::SUCCESS)
    mooseError("Error setting up NEML model history!");
 
  _linear_rotation[_qp].zero();

  _energy[_qp] = 0.0;
  _dissipation[_qp] = 0.0;
}

void tensor_neml(const RankTwoTensor &in, double *const out) {
  double inds[6][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};
  double mults[6] = {1.0, 1.0, 1.0, sqrt(2.0), sqrt(2.0), sqrt(2.0)};

  for (int i = 0; i < 6; i++) {
    out[i] = in(inds[i][0], inds[i][1]) * mults[i];
  }
}

void neml_tensor(const double *const in, RankTwoTensor &out) {
  double inds[6][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};
  double mults[6] = {1.0, 1.0, 1.0, sqrt(2.0), sqrt(2.0), sqrt(2.0)};

  for (int i = 0; i < 6; i++) {
    out(inds[i][0], inds[i][1]) = in[i] / mults[i];
    out(inds[i][1], inds[i][0]) = in[i] / mults[i];
  }
}

void neml_tangent(const double *const in, RankFourTensor &out) {
  double inds[6][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};
  double mults[6] = {1.0, 1.0, 1.0, sqrt(2.0), sqrt(2.0), sqrt(2.0)};

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      out(inds[i][0], inds[i][1], inds[j][0], inds[j][1]) =
          in[i * 6 + j] / (mults[i] * mults[j]);
      out(inds[i][1], inds[i][0], inds[j][0], inds[j][1]) =
          in[i * 6 + j] / (mults[i] * mults[j]);
      out(inds[i][0], inds[i][1], inds[j][1], inds[j][0]) =
          in[i * 6 + j] / (mults[i] * mults[j]);
      out(inds[i][1], inds[i][0], inds[j][1], inds[j][0]) =
          in[i * 6 + j] / (mults[i] * mults[j]);
    }
  }
}

void tensor_skew(const RankTwoTensor &in, double *const out) {
  out[0] = -in(1, 2);
  out[1] = in(0, 2);
  out[2] = -in(0, 1);
}

void skew_tensor(const double *const in, RankTwoTensor &out) {
  out.zero();
  out(0, 1) = -in[2];
  out(0, 2) = in[1];
  out(1, 0) = in[2];
  out(1, 2) = -in[0];
  out(2, 0) = -in[1];
  out(2, 1) = in[0];
}

void recombine_tangent(const double *const Dpart, const double *const Wpart,
                       RankFourTensor &out) {
  std::vector<double> data(81);
  neml::transform_fourth(Dpart, Wpart, &data[0]);
  out.fillFromInputVector(data, RankFourTensor::FillMethod::general);
}
