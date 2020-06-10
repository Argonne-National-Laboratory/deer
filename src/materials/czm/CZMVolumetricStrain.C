#include "CZMVolumetricStrain.h"
#include "RankFourTensor.h"
#include "RotationMatrix.h"

registerMooseObject("DeerApp", CZMVolumetricStrain);

InputParameters CZMVolumetricStrain::validParams() {
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription(
      "Class for computing the interface volume and volume rate");
  params.addRequiredCoupledVar(
      "displacements",
      "The string of displacements suitable for the problem statement");
  params.addClassDescription("Base czm material for large deformation");
  params.addParam<bool>("large_kinematics", true,
                        "Use large displacement kinematics.");
  // params.addParam<bool>("use_area_change", true,
  //                       "account for interface area changes");

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
      _large_kinematics(getParam<bool>("large_kinematics")) {}

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

  RankTwoTensor F_average = (F + F_neighbor) * .5;
  RankTwoTensor F_average_old = (F_old + F_neighbor_old) * .5;

  /// compute needed rotations
  RankTwoTensor R_avg, R_avg_old, DR_avg;
  F_average.getRUDecompositionRotation(R_avg);
  F_average_old.getRUDecompositionRotation(R_avg_old);
  DR_avg = R_avg - R_avg_old;

  RealVectorValue n_average = R_avg * _normals[_qp];
  RealVectorValue Dn_average = DR_avg * _normals[_qp];

  RealVectorValue jump, Djump;

  for (unsigned int i = 0; i < 3; i++) {
    jump(i) = (*_disp_neighbor[i])[_qp] - (*_disp[i])[_qp];
    Djump(i) = jump(i) - ((*_disp_neighbor_old[i])[_qp] - (*_disp_old[i])[_qp]);
  }
  Real dadA = F.det() * (F.inverse().transpose() * _normals[_qp]).norm();
  Real dadA_neighbor =
      F_neighbor.det() *
      (-F_neighbor.inverse().transpose() * _normals[_qp]).norm();
  Real dadA_average = (dadA + dadA_neighbor) * 0.5;

  Real dadA_old =
      F_old.det() * (F_old.inverse().transpose() * _normals[_qp]).norm();
  Real dadA_neighbor_old =
      F_neighbor_old.det() *
      (-F_neighbor_old.inverse().transpose() * _normals[_qp]).norm();
  Real dadA_average_old = (dadA_old + dadA_neighbor_old) * 0.5;
  Real area_increment = dadA_average - dadA_average_old;

  RankTwoTensor uinc_outer_n = outer_product(Djump, n_average);
  RankTwoTensor u_outer_ninc = outer_product(jump, Dn_average);
  RankTwoTensor u_outer_n = outer_product(jump, n_average);
  RankTwoTensor n_outer_n = outer_product(n_average, n_average);
  RankTwoTensor n_outer_Dn = outer_product(n_average, Dn_average);

  /// down here we compute the total volume increment starting from
  /// V = (jump_i*n_j+jump_j*n_i)*A and taking time derivaitives
  _czm_total_strain_rate[_qp] = (uinc_outer_n + uinc_outer_n.transpose());
  _czm_total_strain_rate[_qp] += (u_outer_ninc + u_outer_ninc.transpose());
  _czm_total_strain_rate[_qp] *= dadA_average;
  _czm_total_strain_rate[_qp] +=
      (u_outer_n + u_outer_n.transpose()) * area_increment;
  _czm_total_strain_rate[_qp] *= 0.5 / _dt;

  /// down here we compute the normal volume increment starting from
  /// V = (jump_i*n_i)*(n_i*n_j)*A and taking time derivaitives
  _czm_normal_strain_rate[_qp] =
      (Djump * n_average + jump * Dn_average) * n_outer_n;
  _czm_normal_strain_rate[_qp] +=
      u_outer_n * (n_outer_Dn + n_outer_Dn.transpose());
  _czm_normal_strain_rate[_qp] *= dadA_average;
  _czm_normal_strain_rate[_qp] +=
      (jump * n_average * n_outer_n) * (dadA_average - dadA_average_old);
  _czm_normal_strain_rate[_qp] /= _dt;

  /// compute the sliding contribution as the difference between total and
  /// normal
  _czm_sliding_strain_rate[_qp] =
      _czm_total_strain_rate[_qp] - _czm_normal_strain_rate[_qp];
}
