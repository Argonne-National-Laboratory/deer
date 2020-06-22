#include "ComputeNEMLStrainBase.h"

#include "HomogenizationConstraintKernel.h"

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
      _inv_def_grad(declareProperty<RankTwoTensor>("inv_def_grad")),
      _detJ(declareProperty<Real>("detJ")),
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
  if ((_num_hvars != 0) && (_num_hvars !=
                            HomogenizationConstants::required.at(_ld)[_ndisp-1])) {
    mooseError("Strain calculator must either have 0 or ",
               HomogenizationConstants::required.at(_ld)[_ndisp-1],
               " homogenization scalar variables");
  }

  unsigned int total = 9;
  _homogenization_vals.resize(total);

  unsigned int i;
  for (i = 0; i < _num_hvars; i++) {
    _homogenization_vals[i] = &coupledScalarValue("homogenization_variables", i);
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
  _inv_def_grad[_qp] = RankTwoTensor::Identity();
  _detJ[_qp] = 1.0;
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
    if ((_ndisp == 1) || (_ndisp == 3)) {
      return RankTwoTensor((*_homogenization_vals[0])[0],
                           (*_homogenization_vals[1])[0],
                           (*_homogenization_vals[2])[0],
                           (*_homogenization_vals[3])[0],
                           (*_homogenization_vals[4])[0],
                           (*_homogenization_vals[5])[0],
                           (*_homogenization_vals[6])[0],
                           (*_homogenization_vals[7])[0],
                           (*_homogenization_vals[8])[0]);
    }
    else {
      return RankTwoTensor((*_homogenization_vals[0])[0],
                           (*_homogenization_vals[2])[0],
                           (*_homogenization_vals[4])[0],
                           (*_homogenization_vals[3])[0],
                           (*_homogenization_vals[1])[0],
                           (*_homogenization_vals[5])[0],
                           (*_homogenization_vals[6])[0],
                           (*_homogenization_vals[7])[0],
                           (*_homogenization_vals[8])[0]);
    }
  }
  else {
    if (_ndisp == 1) {
      return RankTwoTensor((*_homogenization_vals[0])[0],
                           (*_homogenization_vals[1])[0],
                           (*_homogenization_vals[2])[0],
                           (*_homogenization_vals[3])[0],
                           (*_homogenization_vals[4])[0],
                           (*_homogenization_vals[5])[0],
                           (*_homogenization_vals[6])[0],
                           (*_homogenization_vals[7])[0],
                           (*_homogenization_vals[8])[0]);
    }
    else if (_ndisp == 2) {
      return RankTwoTensor((*_homogenization_vals[0])[0],  // 0,0
                           (*_homogenization_vals[3])[0],  // 1,0
                           (*_homogenization_vals[4])[0],  // 2,0
                           (*_homogenization_vals[2])[0],  // 0,1
                           (*_homogenization_vals[1])[0],  // 1,1
                           (*_homogenization_vals[5])[0],  // 2,1
                           (*_homogenization_vals[6])[0],  // 0,2
                           (*_homogenization_vals[7])[0],  // 1,2
                           (*_homogenization_vals[8])[0]); // 2,2
    }
    else {
      return RankTwoTensor((*_homogenization_vals[0])[0],   // 0,0
                           (*_homogenization_vals[6])[0],   // 1,0
                           (*_homogenization_vals[7])[0],   // 2,0
                           (*_homogenization_vals[5])[0],   // 0,1
                           (*_homogenization_vals[1])[0],   // 1,1
                           (*_homogenization_vals[8])[0],   // 2,1
                           (*_homogenization_vals[4])[0],   // 0,2
                           (*_homogenization_vals[3])[0],   // 1,2
                           (*_homogenization_vals[2])[0]);  // 2,2
    }
  }
}

void ComputeNEMLStrainBase::computeQpProperties() {
  _homogenization_contribution[_qp] = homogenizationContribution();
  _def_grad[_qp] = (RankTwoTensor::Identity() +
                    RankTwoTensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp],
                                  (*_grad_disp[2])[_qp]))
      + _homogenization_contribution[_qp];
  _inv_def_grad[_qp] = _def_grad[_qp].inverse();
  _detJ[_qp] = _def_grad[_qp].det();
}
