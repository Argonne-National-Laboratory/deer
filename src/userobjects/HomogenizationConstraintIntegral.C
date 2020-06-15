//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HomogenizationConstraintIntegral.h"

#include "Function.h"

registerMooseObject("DeerApp", HomogenizationConstraintIntegral);

InputParameters
HomogenizationConstraintIntegral::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addRequiredCoupledVar("homogenization_variables", "The scalar "
                               "variables with the extra gradient components");
  params.addRequiredParam<std::vector<unsigned int>>("constraint_types",
    "Type of each constraint: strain (0) or stress (1)"); 
  params.addRequiredParam<std::vector<FunctionName>>("targets",
                                        "Functions giving the targets to hit");

  return params;
}

HomogenizationConstraintIntegral::HomogenizationConstraintIntegral(const
                                                                   InputParameters
                                                                   & parameters)
  : ElementUserObject(parameters),
    _num_hvars(coupledScalarComponents("homogenization_variables")),
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
HomogenizationConstraintIntegral::initialize()
{
  _residual.resize(_num_hvars);
  _jacobian.resize(_num_hvars);
  for (_h = 0; _h < _num_hvars; _h++) {
    _residual[_h] = 0.0;
    _jacobian[_h].zero();
  }
}

void
HomogenizationConstraintIntegral::execute()
{
  for (_h = 0; _h < _num_hvars; _h++) {
    for (_qp = 0; _qp < _qrule->n_points(); _qp++) {
      Real dV = _JxW[_qp] * _coord[_qp];
      _residual[_h] += computeResidual() * dV;
      _jacobian[_h] += computeJacobian() * dV;
    }
  }
}

void
HomogenizationConstraintIntegral::threadJoin(const UserObject & y)
{
  const HomogenizationConstraintIntegral & other = 
      static_cast<const HomogenizationConstraintIntegral &>(y);
  for (_h = 0; _h < _num_hvars; _h++) {
    _residual[_h] += other._residual[_h];
    _jacobian[_h] += other._jacobian[_h];
  }
}

void
HomogenizationConstraintIntegral::finalize()
{
  return;
}

std::vector<Real>
HomogenizationConstraintIntegral::getResidual()
{
  for (_h = 0; _h < _num_hvars; _h++) {
    gatherSum(_residual[_h]);
  }
  return _residual;
}

std::vector<RankTwoTensor>
HomogenizationConstraintIntegral::getJacobian()
{
  for (_h = 0; _h < _num_hvars; _h++) {
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        gatherSum(_jacobian[_h](i,j));
      }
    }
  }
  return _jacobian;
}

Real
HomogenizationConstraintIntegral::computeResidual()
{
  if (_ctypes[_h] == ConstraintType::Stress) {
    return _stress[_qp](_pinds[_h].first,_pinds[_h].second) - 
        _targets[_h]->value(_t, _q_point[_qp]);
  }
  else {
    return _F[_qp](_pinds[_h].first,_pinds[_h].second) - 
        (_sfacts[_h] + _targets[_h]->value(_t, _q_point[_qp]));
  }
}

RankTwoTensor
HomogenizationConstraintIntegral::computeJacobian()
{
  RankTwoTensor a;
  a.zero();
  return a;
}
