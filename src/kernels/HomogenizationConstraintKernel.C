//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HomogenizationConstraintKernel.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"
#include "Function.h"

#include "libmesh/quadrature.h"

registerMooseObject("DeerApp", HomogenizationConstraintKernel);

InputParameters
HomogenizationConstraintKernel::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("homogenization_variables", "The scalar "
                               "variables with the extra gradient components");

  params.addRequiredParam<unsigned int>("component",
                                        "Which direction this kernel acts in");
  params.addRequiredCoupledVar("displacements", "The displacement components");
  
  params.addRequiredParam<std::vector<unsigned int>>("constraint_types",
    "Type of each constraint: strain (0) or stress (1)"); 
  params.addRequiredParam<std::vector<FunctionName>>("targets",
                                        "Functions giving the targets to hit");

  return params;
}

HomogenizationConstraintKernel::HomogenizationConstraintKernel(const InputParameters & parameters)
  : Kernel(parameters),
    _ld(getParam<bool>("use_displaced_mesh")),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_nums(_ndisp),
    _disp_vars(_ndisp),
    _grad_disp(_ndisp),
    _num_hvars(coupledScalarComponents("homogenization_variables")),
    _homogenization_nums(_num_hvars),
    _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
    _material_jacobian(
        getMaterialPropertyByName<RankFourTensor>("material_jacobian")),
    _F(getMaterialPropertyByName<RankTwoTensor>("def_grad"))
{
  const std::vector<FunctionName> & names =
      getParam<std::vector<FunctionName>>("targets");

  unsigned int nfns = names.size();
  if (nfns != _num_hvars) {
    mooseError("Number of target functions must match the number of "
               "homogenization variables");
  }
  for (unsigned int i = 0; i < nfns; i++) {
    const Function * const f = &getFunctionByName(names[i]);
    if (!f) mooseError("Function ", names[i], " not found.");
    _targets.push_back(f);
  }

  const std::vector<unsigned int> & types = 
      getParam<std::vector<unsigned int>>("constraint_types");
  if (types.size() != _num_hvars) {
    mooseError("Number of constraint types must match the number of "
               "homogenization variables");
  }
  for (unsigned int i = 0; i < _num_hvars; i++) {
    if (types[i] == 0) {
      _ctypes.push_back(ConstraintType::Strain);
    }
    else if (types[i] == 1) {
      _ctypes.push_back(ConstraintType::Stress);
    }
    else {
      mooseError("Constraint types must be either 0 (strain) or 1 (stress)");
    }
  }
}

void
HomogenizationConstraintKernel::initialSetup() {
  for (unsigned int i = 0; i < _ndisp; i++) {
    _disp_nums[i] = coupled("displacements", i);
    _disp_vars[i] = getVar("displacements", i);
    _grad_disp[i] = &coupledGradient("displacements", i);
  }

  // Do some checking on the number of homogenization variables
  unsigned int needed = (_ld ? _ndisp*_ndisp : (_ndisp*_ndisp+_ndisp)/2);
  if ((_num_hvars != 0) && (_num_hvars != needed)) {
    mooseError("Strain calculator must either have 0 or ", needed, 
               " homogenization scalar variables");
  }

  for (unsigned int i = 0; i < _num_hvars; i++) {
    _homogenization_nums[i] = coupledScalar("homogenization_variables", i);
  }
}

void
HomogenizationConstraintKernel::computeResidual()
{
  // Only do this once...
  if (_component == 0) {
    // But do it for each variable
    for (_h = 0; _h < _num_hvars; _h++) {
      // The contribution for the scalar kernel (size 1)
      DenseVector<Number> & re_scalar =
          _assembly.residualBlock(_homogenization_nums[_h]);
      for (_qp = 0; _qp < _qrule->n_points(); _qp++) {
        Real dV = _JxW[_qp] * _coord[_qp];
        if (_ctypes[_h] == ConstraintType::Stress) {
          re_scalar(0) += (_stress[_qp](_pinds[_h].first,_pinds[_h].second) -
                           _targets[_h]->value(_t, _q_point[_qp])) * dV;
        }
        else {
          re_scalar(0) += (_F[_qp](_pinds[_h].first,_pinds[_h].second) - 
                           (_sfacts[_h] + _targets[_h]->value(_t, _q_point[_qp]))) * dV;
        }
      }
    }
  }
}

void
HomogenizationConstraintKernel::computeOffDiagJacobianScalar(unsigned int jvar)
{
  for (_h = 0; _h < _num_hvars; _h++) {
    if (jvar == _homogenization_nums[_h]) break;
  }
  if (_h == _num_hvars) return; // Not our scalars apparently

  DenseMatrix<Number> & ken = _assembly.jacobianBlock(_var.number(), jvar);
  DenseMatrix<Number> & kne = _assembly.jacobianBlock(jvar, _var.number());
  MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);
  
  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < jv.order(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      {
        Real dV = _JxW[_qp] * _coord[_qp];
        ken(_i, _j) += 0;
        kne(_j, _i) += 0;
      }
}
