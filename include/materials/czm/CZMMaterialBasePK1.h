//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceMaterial.h"

/**
 * This is the base Material class for implementing a traction separation
 * material model including large rotation and area changes. The responsibility
 * of this class is to rotate the displacement jump from global to local
 * coordinate and rotate back traction and traction derivatives. The local
 * coordinate system assumes the following order: opening, tangential1,
 * tangential2. Note that tangential1, tangential2 are arbitrary and therefore
 * the interface assumes an in-plane isotropic behavior. By overriding
 * computeTractionIncrementAndDerivatives in aderived class, different traction
 * separation laws can be implemented. The
 * computeTractionIncrementAndDerivatives method assumes calculations are
 * performed in the local frame. CZM laws should always be implemented in 3D
 * even if they are going to be used in 2D or 1D simulations. This class assumes
 * large rotation, large area changes, and that the traction separation law is
 * only dependent upon the the displacement jump.
 */
class CZMMaterialBasePK1 : public InterfaceMaterial {

public:
  static InputParameters validParams();
  CZMMaterialBasePK1(const InputParameters &parameters);

protected:
  virtual void initialSetup() override;

  /// this method is responsible for computing the traction incremenet  and its
  /// derivatives w.r.t. the interface dispalcement jump increment. See
  /// PureElasticPK1 for an example
  virtual void computeTractionIncrementAndDerivatives() = 0;

  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// normal to the interface
  const MooseArray<Point> &_normals;

  /// number of displacement components
  const unsigned int _ndisp;

  /// the coupled displacement and neighbor displacement values
  ///@{
  std::vector<const VariableValue *> _disp;
  std::vector<const VariableValue *> _disp_neighbor;
  std::vector<MooseVariable *> _disp_vars;
  std::vector<const VariableValue *> _disp_old;
  std::vector<const VariableValue *> _disp_neighbor_old;
  ///@}

  /// the coupled displacement and neighbor displacement gradient
  ///@{
  std::vector<const VariableGradient *> _grad_disp;
  std::vector<const VariableGradient *> _grad_disp_neighbor;
  ///@}

  /// the displacement jump in global and local coordiante
  ///@{
  MaterialProperty<RealVectorValue> &_displacement_jump_global;
  MaterialProperty<RealVectorValue> &_displacement_jump_global_inc;
  MaterialProperty<RealVectorValue> &_displacement_jump_global_old;
  MaterialProperty<RealVectorValue> &_displacement_jump;
  MaterialProperty<RealVectorValue> &_displacement_jump_inc;
  const MaterialProperty<RealVectorValue> &_displacement_jump_old;
  ///@}

  /// the value of the traction in global and local coordinates
  ///@{
  MaterialProperty<RealVectorValue> &_PK1traction;
  MaterialProperty<RealVectorValue> &_PK1traction_inc;
  const MaterialProperty<RealVectorValue> &_PK1traction_old;
  MaterialProperty<RealVectorValue> &_traction;
  MaterialProperty<RealVectorValue> &_traction_inc;
  const MaterialProperty<RealVectorValue> &_traction_old;
  ///@}

  /// the cauchy traction in global coordinates
  MaterialProperty<RealVectorValue> &_traction_deformed;
  /// the PK1 traction in the natural interface coordiante
  MaterialProperty<RealVectorValue> &_PK1traction_natural;

  /// the traction's derivatives wrt the displacement jump in global and
  /// local coordinates
  ///@{
  MaterialProperty<RankTwoTensor> &_dPK1traction_djumpglobal;
  MaterialProperty<RankThreeTensor> &_dPK1traction_dF;
  MaterialProperty<RankTwoTensor> &_dtraction_djump;
  ///@}

  /// the interface deformation gradient and its delta value
  ///@{
  MaterialProperty<RankTwoTensor> &_F_avg;
  const MaterialProperty<RankTwoTensor> &_F_avg_old;
  MaterialProperty<RankTwoTensor> &_DF_avg;
  ///@}

  /// the interface rotation and its delta value
  ///@{
  MaterialProperty<RankTwoTensor> &_R_avg;
  const MaterialProperty<RankTwoTensor> &_R_avg_old;
  MaterialProperty<RankTwoTensor> &_DR_avg;
  ///@}

  // the interface stretch
  MaterialProperty<RankTwoTensor> &_U_avg;
  /// the interface incremental velocity gradient
  MaterialProperty<RankTwoTensor> &_DL_avg;
  /// the interface normal
  MaterialProperty<RealVectorValue> &_n_avg;
  /// the interface area change rate scaled by da
  MaterialProperty<Real> &_dadot_da_avg;
  /// the interface area change
  MaterialProperty<Real> &_da_dA_avg;

  /// the inverse of the deforamtion gradient
  RankTwoTensor _F_avg_inv;
  /// determinant of _F_avg
  Real _J;

  /// coeffcients used to calcualte the PK1 traction increment
  /// traction_inc = _a*((_B+_C)*_traction + D*_traction_inc)
  ///@{
  Real _a;
  RankTwoTensor _B;
  RankTwoTensor _C;
  RankTwoTensor _D;
  ///@}

  /// if true enables large kinematics. Notice that we are not using
  /// _use_displaced_mesh as we always refer to the initial geometry
  const bool _ld;
  /// if true also consider area changes
  const bool _use_area_change;

  /// method used to calculate the deformation gradient, the velocity gradient
  /// and other necessary kinematic quantities
  void computeFandL();

  /// method computing the derative of the interface rotation w.r.t. the
  /// deformation gradient
  RankFourTensor computedRdF(const RankTwoTensor &R,
                             const RankTwoTensor &U) const;

  RankFourTensor computedFinversedF(const RankTwoTensor &F_inv) const {
    return -RikRlj(F_inv, F_inv);
  };

  /// methods used to perform unsual tensor multiplications in indicail notation
  ///@{

  RankFourTensor RikRlj(const RankTwoTensor &R2,
                        const RankTwoTensor &R2_2) const {
    RankFourTensor res;
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int k = 0; k < 3; k++)
          for (unsigned int l = 0; l < 3; l++)
            res(i, j, k, l) = R2(i, k) * R2_2(l, j);

    return res;
  }

  RankThreeTensor RijklVj(const RankFourTensor &R4,
                          const RealVectorValue &V) const {
    RankThreeTensor res;
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++) {
          res(i, k, l) = 0;
          for (unsigned int j = 0; j < 3; j++)
            res(i, k, l) += R4(i, j, k, l) * V(j);
        }
    return res;
  }

  RankThreeTensor RijklVl(const RankFourTensor &R4,
                          const RealVectorValue &V) const {
    RankThreeTensor res;
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++) {
          res(i, k, l) = 0;
          for (unsigned int j = 0; j < 3; j++)
            res(i, k, l) += R4(i, j, k, l) * V(j);
        }
    return res;
  }

  RankThreeTensor RijklVi(const RankFourTensor &R4,
                          const RealVectorValue &V) const {
    RankThreeTensor res;
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++) {
          res(j, k, l) = 0;
          for (unsigned int i = 0; i < 3; i++)
            res(j, k, l) += R4(i, j, k, l) * V(i);
        }
    return res;
  }

  RankThreeTensor RijRjkl(const RankTwoTensor &R2,
                          const RankThreeTensor &R3) const {
    RankThreeTensor res;
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int k = 0; k < 3; k++)
          for (unsigned int l = 0; l < 3; l++)
            res(i, k, l) += R2(i, j) * R3(j, k, l);
    return res;
  }

  RankFourTensor RijklRjm(const RankFourTensor &R4,
                          const RankTwoTensor &R2) const {
    RankFourTensor res;
    res.zero();
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
          for (unsigned int m = 0; m < 3; m++)
            for (unsigned int j = 0; j < 3; j++)
              res(i, m, k, l) += R4(i, j, k, l) * R2(j, m);
    return res;
  }

  RankTwoTensor RijklRij(const RankFourTensor &R4,
                         const RankTwoTensor &R2) const {
    RankTwoTensor res;
    for (unsigned int k = 0; k < 3; k++)
      for (unsigned int l = 0; l < 3; l++)
        for (unsigned int i = 0; i < 3; i++)
          for (unsigned int j = 0; j < 3; j++)
            res(k, l) += R4(i, j, k, l) * R2(i, j);
    return res;
  }

  RankThreeTensor RjkVi(const RankTwoTensor &R2,
                        const RealVectorValue &V) const {
    RankThreeTensor res;

    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int k = 0; k < 3; k++)
          res(i, j, k) = V(i) * R2(j, k);
    return res;
  }
  ///@}
};
