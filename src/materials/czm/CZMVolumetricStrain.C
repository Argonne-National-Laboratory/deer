#include "CZMVolumetricStrain.h"
#include "RankFourTensor.h"
#include "RotationMatrix.h"

registerMooseObject("DeerApp", CZMVolumetricStrain);

InputParameters CZMVolumetricStrain::validParams() {
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription("Material computing volumetric strain resulting "
                             "from interface separation and sliding."
                             "It also provides the normal and sliding strain "
                             "conitrbuting to the total interface strain.");
  params.addRequiredCoupledVar(
      "displacements",
      "The string of displacements suitable for the problem statement");
  params.addParam<bool>("large_kinematics", true,
                        "If true (default) uses large kinematics to "
                        "properly reorient and scale resulting strains");
  params.set<bool>("use_displaced_mesh", false);

  return params;
}

CZMVolumetricStrain::CZMVolumetricStrain(const InputParameters &parameters)
    : InterfaceMaterial(parameters), _normals(_assembly.normals()),
      _ndisp(coupledComponents("displacements")), _disp(3), _disp_neighbor(3),
      _disp_vars(3), _disp_old(3), _disp_neighbor_old(3),
      _czm_total_strain_rate(
          declareProperty<RankTwoTensor>("czm_total_strain_rate")),
      _czm_normal_strain_rate(
          declareProperty<RankTwoTensor>("czm_normal_strain_rate")),
      _czm_sliding_strain_rate(
          declareProperty<RankTwoTensor>("czm_sliding_strain_rate")),
      _czm_total_strain(declareProperty<RankTwoTensor>("czm_total_strain")),
      _czm_normal_strain(declareProperty<RankTwoTensor>("czm_normal_strain")),
      _czm_sliding_strain(declareProperty<RankTwoTensor>("czm_sliding_strain")),
      _czm_area_mp(declareProperty<Real>("czm_area_mp")),
      _ld(_ndisp == 1 ? false : getParam<bool>("large_kinematics")) {

  // Enforce consistency
  if (_ndisp != _mesh.dimension())
    paramError("displacements",
               "Number of displacements must match problem dimension.");

  // enforce  mataterial formulation
  if (getParam<bool>("use_displaced_mesh"))
    mooseError("CZMVolumetricStrain doesn't work with use_displaced_mesh true");
}

void CZMVolumetricStrain::initialSetup() {
  // initializing the displacement vectors
  for (unsigned int i = 0; i < _ndisp; ++i) {
    _disp[i] = &coupledValue("displacements", i);
    _disp_neighbor[i] = &coupledNeighborValue("displacements", i);
    _disp_vars[i] = getVar("displacements", i);
    _disp_old[i] = &coupledValueOld("displacements", i);
    _disp_neighbor_old[i] = &coupledNeighborValueOld("displacements", i);

    _grad_disp.push_back(&coupledGradient("displacements", i));
    _grad_disp_neighbor.push_back(&coupledNeighborGradient("displacements", i));
    _grad_disp_old.push_back(&coupledGradientOld("displacements", i));
    _grad_disp_neighbor_old.push_back(
        &coupledNeighborGradientOld("displacements", i));
  }

  // All others zero (so this will work naturally for 2D and 1D problems)
  for (unsigned int i = _ndisp; i < 3; i++) {
    _disp[i] = &_zero;
    _disp_neighbor[i] = &_zero;
    _disp_old[i] = &_zero;
    _disp_neighbor_old[i] = &_zero;
    _grad_disp.push_back(&_grad_zero);
    _grad_disp_neighbor.push_back(&_grad_zero);
    _grad_disp_old.push_back(&_grad_zero);
    _grad_disp_neighbor_old.push_back(&_grad_zero);
  }
}

void CZMVolumetricStrain::computeQpProperties() {

  // initialize kinematics variable for small deformation
  _DR_avg = RankTwoTensor::Identity();
  _dadA = 1;
  _Da = 0;
  _n_average = _normals[_qp];
  _Dn_average.zero();

  if (_ld) {
    computeFInterface();
    computeRInterface();

    // compute velocity gradient
    _DL = RankTwoTensor::Identity() - _F_average_old * _F_average.inverse();

    computeNormalInterface();
    computeAreaInterface();
  }
  computeJumpInterface();
  computeInterfaceStrainRates();
  _czm_area_mp[_qp] = _JxW[_qp] * _dadA;
}

void CZMVolumetricStrain::computeJumpInterface() {
  // compute displacement _jump increment
  for (unsigned int i = 0; i < 3; i++) {
    _jump(i) = (*_disp_neighbor[i])[_qp] - (*_disp[i])[_qp];
    _Djump(i) =
        _jump(i) - ((*_disp_neighbor_old[i])[_qp] - (*_disp_old[i])[_qp]);
  }
}

void CZMVolumetricStrain::computeInterfaceStrainRates() {
  // precompute some usefull product
  const RankTwoTensor uinc_outer_n = outer_product(_Djump, _n_average);
  const RankTwoTensor u_outer_ninc = outer_product(_jump, _Dn_average);
  const RankTwoTensor u_outer_n = outer_product(_jump, _n_average);
  const RankTwoTensor n_outer_n = outer_product(_n_average, _n_average);
  const RankTwoTensor n_outer_Dn = outer_product(_n_average, _Dn_average);

  /// here we compute the total volume rate starting from
  /// dV_total = 0.5*(jump_i*n_j+jump_j*n_i)*_Da and taking time derivatives
  _czm_total_strain_rate[_qp] = (uinc_outer_n + uinc_outer_n.transpose());
  _czm_total_strain_rate[_qp] += (u_outer_ninc + u_outer_ninc.transpose());
  _czm_total_strain_rate[_qp] *= _dadA;
  _czm_total_strain_rate[_qp] += (u_outer_n + u_outer_n.transpose()) * _Da;
  _czm_total_strain_rate[_qp] *= 0.5 / _dt;

  /// here we compute the normal volume rate starting from
  /// V = (jump_i*n_i)*(n_i*n_j)*A and taking time derivatives
  _czm_normal_strain_rate[_qp] =
      (_Djump * _n_average + _jump * _Dn_average) * n_outer_n;
  _czm_normal_strain_rate[_qp] +=
      _jump * _n_average * (n_outer_Dn + n_outer_Dn.transpose());
  _czm_normal_strain_rate[_qp] *= _dadA;
  _czm_normal_strain_rate[_qp] += (_jump * _n_average * n_outer_n) * _Da;
  _czm_normal_strain_rate[_qp] /= _dt;

  /// compute the sliding contribution as the difference between total and
  /// normal
  _czm_sliding_strain_rate[_qp] =
      _czm_total_strain_rate[_qp] - _czm_normal_strain_rate[_qp];

  _czm_total_strain[_qp] = 0.5 * (u_outer_n + u_outer_n.transpose());
  _czm_normal_strain[_qp] = _jump * _n_average * n_outer_n;
  _czm_sliding_strain[_qp] = _czm_total_strain[_qp] - _czm_normal_strain[_qp];
}

void CZMVolumetricStrain::computeFInterface() {
  ///  intialize deformation gradients
  RankTwoTensor F = (RankTwoTensor::Identity() +
                     RankTwoTensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp],
                                   (*_grad_disp[2])[_qp]));
  RankTwoTensor F_neighbor = (RankTwoTensor::Identity() +
                              RankTwoTensor((*_grad_disp_neighbor[0])[_qp],
                                            (*_grad_disp_neighbor[1])[_qp],
                                            (*_grad_disp_neighbor[2])[_qp]));

  RankTwoTensor F_old =
      (RankTwoTensor::Identity() + RankTwoTensor((*_grad_disp_old[0])[_qp],
                                                 (*_grad_disp_old[1])[_qp],
                                                 (*_grad_disp_old[2])[_qp]));
  RankTwoTensor F_neighbor_old =
      (RankTwoTensor::Identity() +
       RankTwoTensor((*_grad_disp_neighbor_old[0])[_qp],
                     (*_grad_disp_neighbor_old[1])[_qp],
                     (*_grad_disp_neighbor_old[2])[_qp]));

  _F_average = (F + F_neighbor) * .5;
  _F_average_old = (F_old + F_neighbor_old) * .5;
}

void CZMVolumetricStrain::computeRInterface() {
  /// compute needed rotations
  RankTwoTensor R_avg_old, _DR_avg, _DL;
  _F_average.getRUDecompositionRotation(_R_avg);
  _F_average_old.getRUDecompositionRotation(R_avg_old);
  _DR_avg = _R_avg - R_avg_old;
}

void CZMVolumetricStrain::computeNormalInterface() {
  // compute normal
  _n_average = _R_avg * _normals[_qp];
  // compute interface normal change rate *_dt
  _Dn_average = (_n_average * (_DL * _n_average)) * _n_average -
                _DL.transpose() * _n_average;
}

void CZMVolumetricStrain::computeAreaInterface() {
  // compute area change and area change increment
  _dadA = _F_average.det() *
          (_F_average.inverse().transpose() * _normals[_qp]).norm();
  _Da = (_DL.trace() - _n_average * (_DL * _n_average)) * _dadA;
}
