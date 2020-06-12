#include "ComputeNEMLStrainBase.h"

InputParameters ComputeNEMLStrainBase::validParams() {
  InputParameters params = Material::validParams();

  params.addRequiredCoupledVar("displacements", "Displacement variables");
  params.addParam<bool>("large_kinematics", false,
                        "Use large displacement kinematics.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", "List of eigenstrains to account for.");
  params.suppressParameter<bool>("use_displaced_mesh");

  params.addCoupledVar("homogenization_variables", 0.0, 
                       "The scalar variables providing the homogenization "
                       "contributions.");

  return params;
}

ComputeNEMLStrainBase::ComputeNEMLStrainBase(const InputParameters &parameters)
    : DerivativeMaterialInterface<Material>(parameters),
      _ndisp(coupledComponents("displacements")), _disp(3), _grad_disp(3),
      _strain_inc(declareProperty<RankTwoTensor>("strain_inc")),
      _mechanical_strain_inc(
          declareProperty<RankTwoTensor>("mechanical_strain_inc")),
      _vorticity_inc(declareProperty<RankTwoTensor>("vorticity_inc")),
      _def_grad(declareProperty<RankTwoTensor>("def_grad")),
      _def_grad_old(getMaterialPropertyOld<RankTwoTensor>("def_grad")),
      _df(declareProperty<RankTwoTensor>("df")),
      _eigenstrain_names(
          getParam<std::vector<MaterialPropertyName>>("eigenstrain_names")),
      _eigenstrains(_eigenstrain_names.size()),
      _eigenstrains_old(_eigenstrain_names.size()),
      _ld(getParam<bool>("large_kinematics")),
      _num_hvars(coupledScalarComponents("homogenization_variables")),
      _homogenization_contribution(declareProperty<RankTwoTensor>("homogenization_contribution"))
{
  for (unsigned int i = 0; i < _eigenstrain_names.size(); i++) {
    _eigenstrains[i] =
        &getMaterialProperty<RankTwoTensor>(_eigenstrain_names[i]);
    _eigenstrains_old[i] =
        &getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_names[i]);
  }
  
  // Do some checking on the number of homogenization variables
  unsigned int needed = (_ld ? _ndisp*_ndisp : (_ndisp*_ndisp+_ndisp)/2);
  if ((_num_hvars != 0) && (_num_hvars != needed)) {
    mooseError("Strain calculator must either have 0 or ", needed, 
               " homogenization scalar variables");
  }

  unsigned int total = (_ld ? 9 : 6);
  _homogenization_vals.resize(total);

  unsigned int i;
  for (i = 0; i < _num_hvars; i++) {
    _homogenization_vals[i] = &coupledScalarValue("polarization_stress", i);
  }
  for (; i < total; i++) {
    _homogenization_vals[i] = &_zero;
  }
}

void ComputeNEMLStrainBase::initialSetup() {
  // Enforce consistency
  if (_ndisp != _mesh.dimension()) {
    paramError("displacements",
               "Number of displacements must match problem dimension.");
  }

  for (unsigned int i = 0; i < _ndisp; i++) {
    _disp[i] = &coupledValue("displacements", i);
    _grad_disp[i] = &coupledGradient("displacements", i);
  }

  // All others zero (so this will work naturally for plane strain problems)
  for (unsigned int i = _ndisp; i < 3; i++) {
    _disp[i] = &_zero;
    _grad_disp[i] = &_grad_zero;
  }
}

void ComputeNEMLStrainBase::initQpStatefulProperties() {
  _strain_inc[_qp].zero();
  _mechanical_strain_inc[_qp].zero();
  _vorticity_inc[_qp].zero();
  _def_grad[_qp] = RankTwoTensor::Identity();
  _df[_qp] = RankTwoTensor::Identity();
  _homogenization_contribution[_qp].zero();
}

void ComputeNEMLStrainBase::computeProperties() {
  precalculate();
  DerivativeMaterialInterface<Material>::computeProperties();
}

void ComputeNEMLStrainBase::precalculate() {}

RankTwoTensor ComputeNEMLStrainBase::eigenstrainIncrement() {
  RankTwoTensor res;
  res.zero();
  for (unsigned int i = 0; i < _eigenstrain_names.size(); i++) {
    res += (*_eigenstrains[i])[_qp];
    res -= (*_eigenstrains_old[i])[_qp];
  }

  return res;
}

RankTwoTensor 
ComputeNEMLStrainBase::homogenizationContribution()
{
  if (_ld) {
    return RankTwoTensor((*_homogenization_vals[0])[_qp],
                         (*_homogenization_vals[1])[_qp],
                         (*_homogenization_vals[2])[_qp],
                         (*_homogenization_vals[3])[_qp],
                         (*_homogenization_vals[4])[_qp],
                         (*_homogenization_vals[5])[_qp],
                         (*_homogenization_vals[6])[_qp],
                         (*_homogenization_vals[7])[_qp],
                         (*_homogenization_vals[8])[_qp]);
  }
  else {
    return RankTwoTensor((*_homogenization_vals[0])[_qp],
                         (*_homogenization_vals[1])[_qp],
                         (*_homogenization_vals[2])[_qp],
                         (*_homogenization_vals[3])[_qp],
                         (*_homogenization_vals[4])[_qp],
                         (*_homogenization_vals[5])[_qp]);
  }
}

void ComputeNEMLStrainBase::computeQpProperties() {
  _homogenization_contribution[_qp] = homogenizationContribution();
  _def_grad[_qp] = (RankTwoTensor::Identity() +
                    RankTwoTensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp],
                                  (*_grad_disp[2])[_qp]))
      + _homogenization_contribution[_qp];
}
