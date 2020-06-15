//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HomogenizationConstraintScalarKernel.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableScalar.h"
#include "Function.h"

registerMooseObject("DeerApp", HomogenizationConstraintScalarKernel);

InputParameters
HomogenizationConstraintScalarKernel::validParams()
{
  InputParameters params = ScalarKernel::validParams();
  
  params.addRequiredParam<unsigned int>("ndim", "Number of problem dimensions");
  params.addRequiredCoupledVar("homogenization_variables", "The scalar "
                               "variables with the extra gradient components");
  params.addRequiredParam<unsigned int>("component", 
                                        "The # of the constraint variable");
  params.addRequiredParam<UserObjectName>("integrator",
                                          "The integrator user object doing the "
                                          "element calculations.");

  return params;
}

HomogenizationConstraintScalarKernel::HomogenizationConstraintScalarKernel(const InputParameters & parameters)
  : ScalarKernel(parameters),
    _ndisp(getParam<unsigned int>("ndim")),
    _num_hvars(coupledScalarComponents("homogenization_variables")),
    _homogenization_nums(_num_hvars),
    _h(getParam<unsigned int>("component")),
    _integrator(getUserObject<HomogenizationConstraintIntegral>("integrator"))
{
  for (unsigned int i = 0; i < _num_hvars; i++) {
    _homogenization_nums[i] = coupledScalar("homogenization_variables", i);
  }
  // A basic sanity check
  if ((_h < 0) || (_h >= _num_hvars)) {
    mooseError("The homogenization variable number must be in between ",
               0, " and ", _num_hvars);
  }

  _pinds = _bpinds[_ndisp-1];
}

void
HomogenizationConstraintScalarKernel::reinit()
{
}

void
HomogenizationConstraintScalarKernel::computeResidual()
{
  DenseVector<Number> & re =  _assembly.residualBlock(_var.number());
  Real value = _integrator.getResidual(_h);
  for (_i = 0; _i < re.size(); _i++)
    re(_i) += value; 
}

void
HomogenizationConstraintScalarKernel::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  RankTwoTensor value = _integrator.getJacobian(_h);
  ke(0, 0) += value(_pinds[_h].first, _pinds[_h].second);
}

void 
HomogenizationConstraintScalarKernel::computeOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int k = 0; k < _num_hvars; k++) {
    if (_homogenization_nums[k] == jvar) {
      RankTwoTensor value = _integrator.getJacobian(_h);
      DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);
      ke(0,0) += value(_pinds[k].first, _pinds[k].second);
    }
  }
}
