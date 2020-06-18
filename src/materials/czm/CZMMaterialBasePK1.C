//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMMaterialBasePK1.h"
#include "DeformationGradientTools.h"
#include "RankFourTensor.h"
#include "RotationMatrix.h"

InputParameters CZMMaterialBasePK1::validParams() {
  InputParameters params = InterfaceMaterial::validParams();
  params.addRequiredCoupledVar(
      "displacements",
      "The string of displacements suitable for the problem statement");
  params.addClassDescription(
      "Base czm material for incremtnal material models");
  params.addParam<bool>("large_kinematics", true,
                        "Use large displacement kinematics.");
  params.addParam<bool>("use_area_change", true,
                        "Account for interface area changes");

  return params;
}

CZMMaterialBasePK1::CZMMaterialBasePK1(const InputParameters &parameters)
    : InterfaceMaterial(parameters), _normals(_assembly.normals()),
      _ndisp(coupledComponents("displacements")), _disp(3), _disp_neighbor(3),
      _disp_vars(3), _disp_old(3), _disp_neighbor_old(3),
      _displacement_jump_global(
          declareProperty<RealVectorValue>("displacement_jump_global")),
      _displacement_jump_global_inc(
          declareProperty<RealVectorValue>("displacement_jump_global_inc")),
      _displacement_jump_global_old(
          declareProperty<RealVectorValue>("displacement_jump_global_old")),
      _displacement_jump(declareProperty<RealVectorValue>("displacement_jump")),
      _displacement_jump_inc(
          declareProperty<RealVectorValue>("displacement_jump_inc")),
      _displacement_jump_old(
          getMaterialPropertyOld<RealVectorValue>("displacement_jump")),
      _PK1traction(declareProperty<RealVectorValue>("PK1traction")),
      _PK1traction_inc(declareProperty<RealVectorValue>("PK1traction_inc")),
      _PK1traction_old(getMaterialPropertyOld<RealVectorValue>("PK1traction")),
      _traction(declareProperty<RealVectorValue>("traction")),
      _traction_inc(declareProperty<RealVectorValue>("traction_inc")),
      _traction_old(getMaterialPropertyOld<RealVectorValue>("traction")),
      _traction_deformed(declareProperty<RealVectorValue>("traction_deformed")),
      _PK1traction_natural(
          declareProperty<RealVectorValue>("PK1traction_natural")),
      _dPK1traction_djumpglobal(
          declareProperty<RankTwoTensor>("dPK1traction_djumpglobal")),
      _dPK1traction_dF(declareProperty<RankThreeTensor>("dPK1traction_dF")),
      _dtraction_djump(declareProperty<RankTwoTensor>("dtraction_djump")),

      _F_avg(declareProperty<RankTwoTensor>("F_avg")),
      _F_avg_old(getMaterialPropertyOld<RankTwoTensor>("F_avg")),
      _DF_avg(declareProperty<RankTwoTensor>("DF_avg")),

      _R_avg(declareProperty<RankTwoTensor>("R_avg")),
      _R_avg_old(getMaterialPropertyOld<RankTwoTensor>("R_avg")),
      _DR_avg(declareProperty<RankTwoTensor>("DR_avg")),

      _U_avg(declareProperty<RankTwoTensor>("U_avg")),

      _DL_avg(declareProperty<RankTwoTensor>("DL_avg")),
      _n_avg(declareProperty<RealVectorValue>("n_avg")),
      _dadot_da_avg(declareProperty<Real>("dadot_da_avg")),
      _da_dA_avg(declareProperty<Real>("da_dA_avg")),

      _ld(_ndisp > 1 ? getParam<bool>("large_kinematics") : false),
      _use_area_change(_ndisp > 1 ? getParam<bool>("use_area_change") : false) {
  // Enforce consistency
  if (_ndisp != _mesh.dimension())
    paramError("displacements",
               "Number of displacements must match problem dimension.");
  if (!_ld && _use_area_change)
    mooseError("If use_area_change=true then large_kinematics must be true");
  if (_ld && _ndisp < 2)
    mooseError("Large deforamtion can't work in 1D or with a single "
               "displacement variable");
}

void CZMMaterialBasePK1::initialSetup() {
  // initializing the displacement vectors
  for (unsigned int i = 0; i < _ndisp; ++i) {
    _disp[i] = &coupledValue("displacements", i);
    _disp_neighbor[i] = &coupledNeighborValue("displacements", i);
    _disp_vars[i] = getVar("displacements", i);
    _disp_old[i] = &coupledValueOld("displacements", i);
    _disp_neighbor_old[i] = &coupledNeighborValueOld("displacements", i);

    _grad_disp.push_back(&coupledGradient("displacements", i));
    _grad_disp_neighbor.push_back(&coupledNeighborGradient("displacements", i));
  }

  // All others zero (so this will work naturally for 2D and 1D problems)
  for (unsigned int i = _ndisp; i < 3; i++) {
    _disp[i] = &_zero;
    _disp_neighbor[i] = &_zero;
    _disp_old[i] = &_zero;
    _disp_neighbor_old[i] = &_zero;
    _grad_disp.push_back(&_grad_zero);
    _grad_disp_neighbor.push_back(&_grad_zero);
  }
}

void CZMMaterialBasePK1::initQpStatefulProperties() {
  _PK1traction[_qp] = 0;
  _traction[_qp] = 0;
  _displacement_jump[_qp] = 0;
  _F_avg[_qp] = RankTwoTensor::Identity();
  _R_avg[_qp] = RankTwoTensor::Identity();
}

void CZMMaterialBasePK1::computeQpProperties() {

  computeJumpGlobal();
  initKinematicsVariale();
  computeJumpInterface();
  computeTractionIncrementAndDerivatives();
  updateTraction();

  computedTPK1dJumpGlobal();

  if (!_ld)
    _dPK1traction_dF[_qp].zero();
  else
    computedTPK1dF();
}

void CZMMaterialBasePK1::computeF() {
  RankTwoTensor F = (RankTwoTensor::Identity() +
                     RankTwoTensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp],
                                   (*_grad_disp[2])[_qp]));
  RankTwoTensor F_neighbor = (RankTwoTensor::Identity() +
                              RankTwoTensor((*_grad_disp_neighbor[0])[_qp],
                                            (*_grad_disp_neighbor[1])[_qp],
                                            (*_grad_disp_neighbor[2])[_qp]));

  _F_avg[_qp] = 0.5 * (F + F_neighbor);
  _F_avg_inv = _F_avg[_qp].inverse();
}

void CZMMaterialBasePK1::computeRU() {
  _F_avg[_qp].getRUDecompositionRotation(_R_avg[_qp]);
  _U_avg[_qp] = _R_avg[_qp].transpose() * _F_avg[_qp];
  _DF_avg[_qp] = _F_avg[_qp] * _F_avg_old[_qp].inverse();
  _DR_avg[_qp] = _R_avg[_qp] - _R_avg_old[_qp];
}

void CZMMaterialBasePK1::computeLJ() {
  _DL_avg[_qp] = RankTwoTensor::Identity() - _F_avg_old[_qp] * _F_avg_inv;
  _J = _F_avg[_qp].det();
}

void CZMMaterialBasePK1::computeJumpGlobal() {

  // compute the displacement jump increment
  for (unsigned int i = 0; i < 3; i++) {
    _displacement_jump_global[_qp](i) =
        (*_disp_neighbor[i])[_qp] - (*_disp[i])[_qp];
    _displacement_jump_global_old[_qp](i) =
        (*_disp_neighbor_old[i])[_qp] - (*_disp_old[i])[_qp];
  }
  _displacement_jump_global_inc[_qp] =
      _displacement_jump_global[_qp] - _displacement_jump_global_old[_qp];
}
void CZMMaterialBasePK1::initKinematicsVariale() {

  _Q0 = RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));

  // initialize useful quantities
  _n_avg[_qp] = _normals[_qp];
  _dadot_da_avg[_qp] = 0.;
  _da_dA_avg[_qp] = 1.;
  _DR_avg[_qp] = RankTwoTensor::initNone;
  _R_avg[_qp] = RankTwoTensor::Identity();
  _F_avg_inv = RankTwoTensor::Identity();
  _J = 1;

  // initialize large deformation useful terms
  if (_ld) {
    computeF();
    computeRU();
    if (_use_area_change) {
      computeLJ();
      _n_avg[_qp] = _R_avg[_qp] * _normals[_qp];
      _dadot_da_avg[_qp] =
          _DL_avg[_qp].trace() - _n_avg[_qp] * (_DL_avg[_qp] * _n_avg[_qp]);
      _da_dA_avg[_qp] = _J * (_F_avg_inv.transpose() * _normals[_qp]).norm();
    }
  }

  _a = _da_dA_avg[_qp];
  _C = _DR_avg[_qp] * _Q0.transpose();
  _D = _R_avg[_qp] * _Q0.transpose();
  _B = _dadot_da_avg[_qp] * _D;
}

void CZMMaterialBasePK1::computeJumpInterface() {
  // compute local displacement jump
  _displacement_jump_inc[_qp] =
      _C.transpose() * _displacement_jump_global[_qp] +
      _D.transpose() * _displacement_jump_global_inc[_qp];
  _displacement_jump[_qp] =
      _displacement_jump_old[_qp] + _displacement_jump_inc[_qp];
}

void CZMMaterialBasePK1::updateTraction() {

  _traction[_qp] = _traction_old[_qp] + _traction_inc[_qp];
  _traction_deformed[_qp] = _D * _traction[_qp];

  // assemble PK1 traction
  _PK1traction_inc[_qp] = _C * _traction[_qp] + _D * _traction_inc[_qp];
  if (_ld && _use_area_change) {
    _PK1traction_inc[_qp] += _B * _traction[_qp];
    _PK1traction_inc[_qp] *= _a;
  }
  _PK1traction[_qp] = _PK1traction_old[_qp] + _PK1traction_inc[_qp];
  _PK1traction_natural[_qp] = _D.transpose() * _PK1traction[_qp];
}

void CZMMaterialBasePK1::computedTPK1dJumpGlobal() {
  // compute the PK1 traction derivatives w.r.t the displacment jump in global
  // coordinates
  const RankTwoTensor djump_djumpglobal = (_C + _D).transpose();
  const RankTwoTensor dtraction_djumpglobal =
      _dtraction_djump[_qp] * djump_djumpglobal;
  _dPK1traction_djumpglobal[_qp] =
      _a * ((_B + _C + _D) * dtraction_djumpglobal);
}

void CZMMaterialBasePK1::computedTPK1dF() {

  _dR_dF = DeformationGradientTools::computedRdF(_R_avg[_qp], _U_avg[_qp]);
  computedCoefficientsdF();
  assembledTPK1dF();
};

void CZMMaterialBasePK1::computedCdF() {
  _dC_dF = RijklRjm(_dR_dF, _Q0.transpose());
}

void CZMMaterialBasePK1::computedadF(const Real &Fitr_N_norm,
                                     const RealVectorValue &Fitr_N,
                                     const RankTwoTensor &F_itr,
                                     const RankFourTensor &dFinv_dF) {
  RankTwoTensor R2temp;
  for (unsigned int l = 0; l < 3; l++)
    for (unsigned int m = 0; m < 3; m++) {
      R2temp(l, m) = 0;
      for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
          R2temp(l, m) += Fitr_N(i) * dFinv_dF(j, i, l, m) * _normals[_qp](j);

      R2temp(l, m) *= _J / Fitr_N_norm;
    }
  _da_dF = _J * F_itr * Fitr_N_norm + R2temp;
};

void CZMMaterialBasePK1::computedBdF(const Real &Fitr_N_norm,
                                     const RealVectorValue &Fitr_N,
                                     const RankTwoTensor &F_itr,
                                     const RankFourTensor &dFinv_dF) {

  RankTwoTensor R2temp;

  // compute dBdF
  // --- compute the derivatitve of dtrace(L)dF
  // ------start with dDL_dF
  RankFourTensor dDL_dF;
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int l = 0; l < 3; l++)
        for (unsigned int m = 0; m < 3; m++) {
          dDL_dF(i, j, l, m) = 0;
          for (unsigned int k = 0; k < 3; k++)
            dDL_dF(i, j, l, m) -= _F_avg_old[_qp](i, k) * dFinv_dF(k, j, l, m);
        }

  const RankTwoTensor dtraceDL_dF =
      RankTwoTensor::Identity().initialContraction(dDL_dF);

  // compute the derivatitve of the the second part of dadot_da_avg
  // (n_i*DL_ij*n_j)
  for (unsigned int p = 0; p < 3; p++)
    for (unsigned int q = 0; q < 3; q++) {
      R2temp(p, q) = 0;
      for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
          for (unsigned int k = 0; k < 3; k++)
            for (unsigned int l = 0; l < 3; l++)
              R2temp(p, q) +=
                  _dR_dF(i, k, p, q) * _normals[_qp](k) * _DL_avg[_qp](i, j) *
                      _R_avg[_qp](j, l) * _normals[_qp](l) +
                  _dR_dF(j, l, p, q) * _normals[_qp](l) * _DL_avg[_qp](i, j) *
                      _R_avg[_qp](i, k) * _normals[_qp](k) +
                  dDL_dF(i, j, p, q) * _R_avg[_qp](i, k) * _normals[_qp](k) *
                      _R_avg[_qp](j, l) * _normals[_qp](l);
    }

  // assemble ddsdotdS_dF*_D
  _dB_dF = _D.outerProduct(dtraceDL_dF - R2temp) + _dadot_da_avg[_qp] * _dC_dF;
};

void CZMMaterialBasePK1::computedCoefficientsdF() {
  // compute dCdF
  computedCdF();

  _da_dF.zero();
  _dB_dF.zero();

  if (_use_area_change) {
    // dadF
    const RankTwoTensor F_itr = _F_avg_inv.transpose();
    const RealVectorValue Fitr_N = F_itr * _normals[_qp];
    const Real Fitr_N_norm = Fitr_N.norm();
    const RankFourTensor dFinv_dF =
        DeformationGradientTools::computedFinversedF(_F_avg_inv);

    computedadF(Fitr_N_norm, Fitr_N, F_itr, dFinv_dF);
    computedBdF(Fitr_N_norm, Fitr_N, F_itr, dFinv_dF);

  } // end of area changes derivatives}
}

void CZMMaterialBasePK1::assembledTPK1dF() {
  RealVectorValue Vtemp;
  RankTwoTensor R2temp;
  RankThreeTensor R3temp;
  // assemble _dPK1traction_dF as the sum of three contributions: the 1st
  // terms  includes dadA_dF, the 2nd term includes the contributions of
  // dB_dF, dC_dF and dD_dF, the 3rd term includes the contritbution of
  // djump_dF
  //  T1;
  if (_use_area_change) {
    Vtemp = (_B + _C) * _traction[_qp] + _D * _traction_inc[_qp];
    _dPK1traction_dF[_qp] = RjkVi(_da_dF, Vtemp);
  } else
    _dPK1traction_dF[_qp].zero();

  //  T2;
  Vtemp = _traction[_qp] + _traction_inc[_qp];
  RankThreeTensor dT_dF_temp = _a * RijklVj(_dC_dF, Vtemp);
  if (_use_area_change)
    dT_dF_temp += _a * RijklVj(_dB_dF, _traction[_qp]);
  _dPK1traction_dF[_qp] += dT_dF_temp;

  //  T3;
  R2temp = (_B + _C + _D) * _dtraction_djump[_qp];
  Vtemp = _displacement_jump_global[_qp] + _displacement_jump_global_inc[_qp];
  R3temp = RijklVi(_dC_dF, Vtemp);
  _dPK1traction_dF[_qp] += _a * RijRjkl(R2temp, R3temp);
}
