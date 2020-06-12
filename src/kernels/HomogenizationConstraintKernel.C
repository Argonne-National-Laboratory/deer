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
  params.addRequiredCoupledVar("polarization_stress", "The scalar variables");
  params.addRequiredParam<unsigned int>("component",
                                        "Which direction this kernel acts in");
  params.addRequiredCoupledVar("displacements", "The displacement components");
  
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
    _num_polarization(coupledScalarComponents("polarization_stress")),
    _polarization_nums(_num_polarization),
    _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
    _material_jacobian(
        getMaterialPropertyByName<RankFourTensor>("material_jacobian")),
    _df(getMaterialPropertyByName<RankTwoTensor>("df")),
    _base_stress(getMaterialPropertyByName<RankTwoTensor>("base_stress"))
{
  const std::vector<FunctionName> & names =
      getParam<std::vector<FunctionName>>("targets");
  unsigned int nfns = names.size();
  if (nfns != _num_polarization) {
    mooseError("Number of target functions must match the number of "
               "polarization stress components");
  }
  for (unsigned int i = 0; i < nfns; i++) {
    const Function * const f = &getFunctionByName(names[i]);
    if (!f) mooseError("Function ", names[i], " not found.");
    _targets.push_back(f);
  }
}

void
HomogenizationConstraintKernel::initialSetup() {
  for (unsigned int i = 0; i < _ndisp; i++) {
    _disp_nums[i] = coupled("displacements", i);
    _disp_vars[i] = getVar("displacements", i);
    _grad_disp[i] = &coupledGradient("displacements", i);
  }

  // Do some checking on the provided number of polarization variables
  if ((_num_polarization != 0) && 
      (_num_polarization != (_ndisp*_ndisp+_ndisp)/2)) {
    mooseError("Number of polarization variables should be either zero or "
               "enough to fill in a symmetric nxn tensor (1, 3, 6)");
  }

  for (unsigned int i = 0; i < _num_polarization; i++) {
    _polarization_nums[i] = coupledScalar("polarization_stress", i);
  }
}

void
HomogenizationConstraintKernel::computeResidual()
{
  // Only do this once...
  if (_component == 0) {
    // But do it for each polarization stress
    for (_ps = 0; _ps < _num_polarization; _ps++) {
      // The contribution for the scalar kernel (size 1)
      DenseVector<Number> & re_scalar =
          _assembly.residualBlock(_polarization_nums[_ps]);
      for (_qp = 0; _qp < _qrule->n_points(); _qp++) {
        Real dV = _JxW[_qp] * _coord[_qp];
        re_scalar(0) += (_stress[_qp](_pinds[_ps].first,_pinds[_ps].second) -
                         _targets[_ps]->value(_t, _q_point[_qp])) * dV;
      }
    }
  }
}

void
HomogenizationConstraintKernel::computeOffDiagJacobianScalar(unsigned int jvar)
{
  for (_ps = 0; _ps < _num_polarization; _ps++) {
    if (jvar == _polarization_nums[_ps]) break;
  }
  if (_ps == _num_polarization) return; // Not our scalars apparently

  DenseMatrix<Number> & ken = _assembly.jacobianBlock(_var.number(), jvar);
  DenseMatrix<Number> & kne = _assembly.jacobianBlock(jvar, _var.number());
  MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);
  
  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < jv.order(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      {
        Real dV = _JxW[_qp] * _coord[_qp];
        ken(_i, _j) += computeQpDispJacobian() * dV;
        kne(_j, _i) += computeQpScalarJacobian() * dV;
      }
}

Real
HomogenizationConstraintKernel::computeQpDispJacobian()
{
  if (_component == _pinds[_ps].first) {
    return _grad_test[_i][_qp](_pinds[_ps].second);
  }
  else {
    return 0;
  }
}

Real
HomogenizationConstraintKernel::computeQpScalarJacobian()
{
  // There some funny "we are a galerkin method" going on here
  size_t i = _pinds[_ps].first;
  size_t j = _pinds[_ps].second;
  size_t k = _component;
  Real value = 0.0;
  for (size_t l = 0; l < _ndisp; l++) {
    value += _material_jacobian[_qp](i,j,k,l) * _grad_phi[_i][_qp](l);
  }

  return 0.0;
}
