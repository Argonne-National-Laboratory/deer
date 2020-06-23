#include "TotalStressDivergenceNEML.h"

#include "HomogenizationConstraintIntegral.h"

registerMooseObject("DeerApp", TotalStressDivergenceNEML);

InputParameters TotalStressDivergenceNEML::validParams() {
  InputParameters params = Kernel::validParams();

  params.addRequiredParam<unsigned int>("component",
                                        "Which direction this kernel acts in");
  params.addRequiredCoupledVar("displacements", "The displacement components");
  params.addParam<bool>("large_kinematics", false,
                        "Use large displacement kinematics.");
  params.suppressParameter<bool>("use_displaced_mesh");

  params.addCoupledVar("macro_gradient",
                       "Optional scalar field with the macro gradient");
  params.addParam<std::vector<unsigned int>>("constraint_types",
    "Type of each constraint: strain (0) or stress (1)");

  return params;
}

TotalStressDivergenceNEML::TotalStressDivergenceNEML(const InputParameters &parameters)
    : DerivativeMaterialInterface<Kernel>(parameters),
      _ld(getParam<bool>("large_kinematics")),
      _component(getParam<unsigned int>("component")),
      _ndisp(coupledComponents("displacements")), _disp_nums(_ndisp),
      _disp_vars(_ndisp), _grad_disp(_ndisp),
      _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
      _material_jacobian(
          getMaterialPropertyByName<RankFourTensor>("material_jacobian")),
      _inv_def_grad(getMaterialPropertyByName<RankTwoTensor>("inv_def_grad")),
      _detJ(getMaterialPropertyByName<Real>("detJ")),
      _df(getMaterialPropertyByName<RankTwoTensor>("df")),
      _macro_gradient_num(isCoupledScalar("macro_gradient", 0) ?
                          coupledScalar("macro_gradient") : 0),
      _macro_gradient(isCoupledScalar("macro_gradient", 0) ?
                      getScalarVar("macro_gradient", 0) : nullptr),
      _indices(HomogenizationConstants::indices.at(_ld)[_ndisp-1])
{
  if (isCoupledScalar("macro_gradient", 0)) {
    // Check the order of the scalar variable
    unsigned int needed = HomogenizationConstants::required.at(_ld)[_ndisp-1];
    if (_macro_gradient->order() != needed)
      mooseError("The homogenization macro gradient variable must have order ",
                 needed);
    
    // Check the number of constraints
    auto types = getParam<std::vector<unsigned int>>("constraint_types");
    if (types.size() != needed)
      mooseError("The kernel must be supplied ", needed, " constraint types.");
  }
}

void 
TotalStressDivergenceNEML::initialSetup() {
  for (unsigned int i = 0; i < _ndisp; i++) {
    _disp_nums[i] = coupled("displacements", i);
    _disp_vars[i] = getVar("displacements", i);
    _grad_disp[i] = &coupledGradient("displacements", i);
  }
  
  if (isCoupledScalar("macro_gradient", 0)) {
    auto types = getParam<std::vector<unsigned int>>("constraint_types");
    for (unsigned int i = 0; i < types.size(); i++) {
      if (types[i] == 0)
        _ctypes.push_back(ConstraintType::Strain);
      else if (types[i] == 1)
        _ctypes.push_back(ConstraintType::Stress);
      else
        mooseError("Constraint types need to be 0 or 1");
    }
  }
}

Real TotalStressDivergenceNEML::computeQpResidual() {
  if (_ld) 
    return largeDeformationResidual(_grad_test[_i][_qp]);
  else 
    return smallDeformationResidual(_grad_test[_i][_qp]);
}

Real
TotalStressDivergenceNEML::largeDeformationResidual(const RealGradient & grad_phi)
{
  return _detJ[_qp] * _stress[_qp].row(_component) * 
      (_inv_def_grad[_qp].transpose() * grad_phi);
}

Real
TotalStressDivergenceNEML::smallDeformationResidual(const RealGradient & grad_phi)
{
  return _stress[_qp].row(_component) * grad_phi;
}

Real TotalStressDivergenceNEML::computeQpJacobian() 
{
  Real value = 0.0;
  
  if (_ld) {
    value += largeDeformationMatJac(_component, _component,
                                    _grad_test[_i][_qp],
                                    _grad_phi[_j][_qp]);
    value += largeDeformationGeoJac(_component, _component,
                                    _grad_test[_i][_qp],
                                    _grad_phi[_j][_qp]);
  }
  else {
    value += smallDeformationMatJac(_component, _component,
                                    _grad_test[_i][_qp],
                                    _grad_phi[_j][_qp]);
  }

  return value;
}

Real TotalStressDivergenceNEML::computeQpOffDiagJacobian(unsigned int jvar) {
  Real value = 0.0;

  for (unsigned int cc = 0; cc < _ndisp; cc++) {
    if (jvar == _disp_nums[cc]) {
      if (_ld) {
        value += largeDeformationMatJac(_component, cc,
                                        _grad_test[_i][_qp],
                                        _disp_vars[cc]->gradPhi()[_j][_qp]);
        value += largeDeformationGeoJac(_component, cc,
                                        _grad_test[_i][_qp],
                                        _disp_vars[cc]->gradPhi()[_j][_qp]);
      }
      else {
        value += smallDeformationMatJac(_component, cc,
                                        _grad_test[_i][_qp],
                                        _disp_vars[cc]->gradPhi()[_j][_qp]);
      }
      break;
    }
  }

  return value;
}

Real 
TotalStressDivergenceNEML::smallDeformationMatJac(
    unsigned int i, unsigned int k, const RealGradient & grad_phi,
    const RealGradient & grad_psi)
{
  Real value = 0.0;
  for (unsigned int j = 0; j < _ndisp; j++) 
    for (unsigned int l = 0; l < _ndisp; l++) 
      value += _material_jacobian[_qp](i,j,k,l) * grad_phi(j) * grad_psi(l);
    
  return value;
}

Real
TotalStressDivergenceNEML::largeDeformationMatJac(
    unsigned int i, unsigned int k, const RealGradient & grad_phi,
    const RealGradient & grad_psi)
{
  auto GPhi = _inv_def_grad[_qp].transpose() * grad_phi;
  auto GPsi = _inv_def_grad[_qp].transpose() * grad_psi;

  Real value = 0.0;

  for (unsigned int j = 0; j < _ndisp; j++) 
    for (unsigned int m = 0; m < _ndisp; m++) 
      for (unsigned int n = 0; n < _ndisp; n++) 
        value += _detJ[_qp] * _material_jacobian[_qp](i,j,m,n) * 
            _df[_qp](m,k) * GPsi(n) * GPhi(j);
      
  return value;
}

Real
TotalStressDivergenceNEML::largeDeformationGeoJac(
    unsigned int i, unsigned int k, const RealGradient & grad_phi,
    const RealGradient & grad_psi)
{
  auto GPhi = _inv_def_grad[_qp].transpose() * grad_phi;
  auto GPsi = _inv_def_grad[_qp].transpose() * grad_psi;
  
  Real value = 0.0;

  for (unsigned int j = 0; j < _ndisp; j++) 
    value += _detJ[_qp] * _stress[_qp](i,j) * 
        (GPsi(k) * GPhi(j) - GPsi(j) * GPhi(k));

  return value;
}

void
TotalStressDivergenceNEML::computeOffDiagJacobianScalar(unsigned int jvar)
{
  if (jvar == _macro_gradient_num) {
    DenseMatrix<Number> & ken = _assembly.jacobianBlock(_var.number(), jvar);
    DenseMatrix<Number> & kne = _assembly.jacobianBlock(jvar, _var.number()); // comment this out for FDP checking
    
    for (_qp = 0; _qp < _qrule->n_points(); _qp++) {
      Real dV = _JxW[_qp] * _coord[_qp];
      for (_h = 0; _h < _macro_gradient->order(); _h++) {
        for (_i = 0; _i < _test.size(); _i++) {
          ken(_i, _h) += computeBaseJacobian() * dV;
          kne(_h, _i) += computeConstraintJacobian() * dV; // comment this out for FDP checking
        }
      }
    }
  }
}

Real
TotalStressDivergenceNEML::computeBaseJacobian()
{ 
  Real value = 0.0;
  for (unsigned int j = 0; j < _ndisp; j++) 
    value += _material_jacobian[_qp](_component, j, _indices[_h].first,
                                     _indices[_h].second) *
        _grad_test[_i][_qp](j);

  return value;
}

Real
TotalStressDivergenceNEML::computeConstraintJacobian()
{
  Real value = 0.0;
  if (_ctypes[_h] == ConstraintType::Stress) {
    // Seems right
    for (unsigned int l = 0; l < _ndisp; l++) {
      value += _material_jacobian[_qp](_indices[_h].first, _indices[_h].second,
                                       _component, l) * _grad_phi[_i][_qp](l);
    }
  }
  else {
    for (unsigned int l = 0; l < _ndisp; l++) {
      if ((_indices[_h].first == _component) && (_indices[_h].second == l))
        value += 0.5*_grad_phi[_i][_qp](l);
      if ((_indices[_h].second == _component) && (_indices[_h].first == l))
        value += 0.5*_grad_phi[_i][_qp](l);
    }
  }
  return value;
}
