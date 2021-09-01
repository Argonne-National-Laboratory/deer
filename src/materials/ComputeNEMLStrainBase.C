#include "ComputeNEMLStrainBase.h"

#include "HomogenizationConstraintIntegral.h" // Index data

InputParameters ComputeNEMLStrainBase::validParams() {
  InputParameters params = Material::validParams();

  params.addRequiredCoupledVar("displacements", "Displacement variables");
  params.addParam<bool>("large_kinematics", false,
                        "Use large displacement kinematics.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", "List of eigenstrains to account for.");
  params.suppressParameter<bool>("use_displaced_mesh");

  params.addCoupledVar("macro_gradient",
                       "Optional scalar field with the macro gradient");

  return params;
}

ComputeNEMLStrainBase::ComputeNEMLStrainBase(const InputParameters &parameters)
    : DerivativeMaterialInterface<Material>(parameters),
      _ndisp(coupledComponents("displacements")), _disp(3), _grad_disp(3),
      _strain_inc(declareProperty<RankTwoTensor>("strain_inc")),
      _mechanical_strain_unrotated_inc(
          declareProperty<RankTwoTensor>("mechanical_strain_unrotated_inc")),
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
      _macro_gradient(isCoupledScalar("macro_gradient", 0) ?
                      coupledScalarValue("macro_gradient") : _zero),
      _homogenization_contribution(declareProperty<RankTwoTensor>("homogenization_contribution"))
{
  for (unsigned int i = 0; i < _eigenstrain_names.size(); i++) {
    _eigenstrains[i] =
        &getMaterialProperty<RankTwoTensor>(_eigenstrain_names[i]);
    _eigenstrains_old[i] =
        &getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_names[i]);
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
  _mechanical_strain_unrotated_inc[_qp].zero();
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
  if (isCoupledScalar("macro_gradient",0)) {
    if (_ld) {
      if (_ndisp == 1) {
        return RankTwoTensor(_macro_gradient[0],
                             0,
                             0,
                             0,
                             0,
                             0,
                             0,
                             0,
                             0);
      }
      else if (_ndisp == 2) {
        return RankTwoTensor(_macro_gradient[0],
                             _macro_gradient[2],
                             0,
                             _macro_gradient[3],
                             _macro_gradient[1],
                             0,
                             0,
                             0,
                             0);
      }
      else {
        return RankTwoTensor(_macro_gradient[0],
                             _macro_gradient[1],
                             _macro_gradient[2],
                             _macro_gradient[3],
                             _macro_gradient[4],
                             _macro_gradient[5],
                             _macro_gradient[6],
                             _macro_gradient[7],
                             _macro_gradient[8]);
      }
    }
    else {
      if (_ndisp == 1) {
        return RankTwoTensor(_macro_gradient[0],
                             0,
                             0,
                             0,
                             0,
                             0,
                             0,
                             0,
                             0);
      }
      else if (_ndisp == 2) {
        return RankTwoTensor(_macro_gradient[0],
                             0,
                             0,
                             _macro_gradient[2],
                             _macro_gradient[1],
                             0,
                             0,
                             0,
                             0);
      }
      else {
        return RankTwoTensor(_macro_gradient[0],
                             0.0,
                             0.0,
                             _macro_gradient[5],
                             _macro_gradient[1],
                             0.0,
                             _macro_gradient[4],
                             _macro_gradient[3],
                             _macro_gradient[2]);
      }
    }
  }
  else {
    RankTwoTensor res;
    res.zero();
    return res;
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
