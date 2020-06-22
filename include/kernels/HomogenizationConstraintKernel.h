//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

namespace HomogenizationConstants
{
  typedef std::vector<std::pair<unsigned int, unsigned int>> index_list;
  const std::map<bool, std::vector<index_list>> indices {
  {true, {
    {{0,0}},
    {{0,0},{1,1},{1,0},{0,1}},
    {{0,0},{1,0},{2,0},{0,1},{1,1},{2,1},{2,0},{2,1},{2,2}}
         }},
  {false, {
    {{0,0}},
    {{0,0},{1,1},{0,1}},
    {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}}
          }}};
}

class HomogenizationConstraintKernel : public Kernel
{
 public:
  static InputParameters validParams();

  HomogenizationConstraintKernel(const InputParameters & parameters);
 
  virtual void initialSetup();
  virtual void computeOffDiagJacobianScalar(unsigned int jvar);

 protected:
  virtual Real computeQpResidual() {return 0;};
  virtual Real computeDisplacementJacobian();
  virtual Real computeConstraintJacobian();

 protected:
  unsigned int _h;
  unsigned int _ii;

  const bool _ld;

  unsigned int _ndisp;

  std::vector<unsigned int> _disp_nums;
  std::vector<MooseVariable *> _disp_vars;
  std::vector<const VariableGradient *> _grad_disp;

  unsigned int _num_hvars;
  std::vector<unsigned int> _homogenization_nums;

  const MaterialProperty<RankTwoTensor> &_stress;
  const MaterialProperty<RankFourTensor> &_material_jacobian;
  const MaterialProperty<RankTwoTensor> &_F;

  const HomogenizationConstants::index_list _indices;

  enum class ConstraintType { Stress, Strain };
  std::vector<ConstraintType> _ctypes;

  std::vector<const Function*> _targets;
};
