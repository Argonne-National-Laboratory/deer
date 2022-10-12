#include "ComputeAnisotropicThermalExpansionEigenstrain.h"

#include "RankTwoTensor.h"

registerMooseObject("DeerApp", ComputeAnisotropicThermalExpansionEigenstrain);

InputParameters
ComputeAnisotropicThermalExpansionEigenstrain::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute the anisotropic thermal expansion eigenstrain given three "
                             "CTE functions and three bases.");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addRequiredParam<MaterialPropertyName>(
      "eigenstrain_name",
      "Material property name for the eigenstrain tensor computed "
      "by this model. IMPORTANT: The name of this property must "
      "also be provided to the strain calculator.");
  params.addRequiredCoupledVar("temperature", "Coupled temperature");
  params.addRequiredCoupledVar("stress_free_temperature",
                               "Reference temperature at which there is no "
                               "thermal expansion for thermal eigenstrain "
                               "calculation");
  params.addRequiredParam<MaterialPropertyName>("basis_0", "The first basis");
  params.addRequiredParam<MaterialPropertyName>("basis_1", "The second basis");
  params.addRequiredParam<MaterialPropertyName>("basis_2", "The third basis");
  params.addRequiredParam<FunctionName>("CTE_0",
                                        "The instantaneous CTE function for the first basis");
  params.addRequiredParam<FunctionName>("CTE_1",
                                        "The instantaneous CTE function for the second basis");
  params.addRequiredParam<FunctionName>("CTE_2",
                                        "The instantaneous CTE function for the third basis");
  return params;
}

ComputeAnisotropicThermalExpansionEigenstrain::ComputeAnisotropicThermalExpansionEigenstrain(
    const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _eigenstrain(declareProperty<RankTwoTensor>(
        _base_name + getParam<MaterialPropertyName>("eigenstrain_name"))),
    _T(coupledValue("temperature")),
    _T_old(coupledValueOld("temperature")),
    _n0(getMaterialProperty<RealVectorValue>(_base_name +
                                             getParam<MaterialPropertyName>("basis_0"))),
    _n1(getMaterialProperty<RealVectorValue>(_base_name +
                                             getParam<MaterialPropertyName>("basis_1"))),
    _n2(getMaterialProperty<RealVectorValue>(_base_name +
                                             getParam<MaterialPropertyName>("basis_2"))),
    _CTE0(getFunction("CTE_0")),
    _CTE1(getFunction("CTE_1")),
    _CTE2(getFunction("CTE_2")),
    _thermal_strain_0(declareProperty<Real>(_base_name + "thermal_strian_0")),
    _thermal_strain_1(declareProperty<Real>(_base_name + "thermal_strian_1")),
    _thermal_strain_2(declareProperty<Real>(_base_name + "thermal_strian_2")),
    _thermal_strain_0_old(getMaterialPropertyOld<Real>(_base_name + "thermal_strian_0")),
    _thermal_strain_1_old(getMaterialPropertyOld<Real>(_base_name + "thermal_strian_1")),
    _thermal_strain_2_old(getMaterialPropertyOld<Real>(_base_name + "thermal_strian_2")),
    _step_zero(declareRestartableData<bool>("step_zero", true))
{
}

void
ComputeAnisotropicThermalExpansionEigenstrain::initQpStatefulProperties()
{
  _eigenstrain[_qp].zero();
  _thermal_strain_0[_qp] = 0;
  _thermal_strain_1[_qp] = 0;
  _thermal_strain_2[_qp] = 0;
}

void
ComputeAnisotropicThermalExpansionEigenstrain::computeQpProperties()
{
  if (_t_step > 0)
    _step_zero = false;

  // Skip the eigenstrain calculation in step zero because no solution is computed during
  // the zeroth step, hence computing the eigenstrain in the zeroth step would result in
  // an incorrect calculation of mechanical_strain, which is stateful.
  if (_step_zero)
    return;

  const Real alpha0 = _CTE0.value(_T[_qp]);
  const Real alpha1 = _CTE1.value(_T[_qp]);
  const Real alpha2 = _CTE2.value(_T[_qp]);

  const Real alpha0_old = _CTE0.value(_T_old[_qp]);
  const Real alpha1_old = _CTE1.value(_T_old[_qp]);
  const Real alpha2_old = _CTE2.value(_T_old[_qp]);

  _thermal_strain_0[_qp] =
      _thermal_strain_0_old[_qp] + (_T[_qp] - _T_old[_qp]) * (alpha0 + alpha0_old) / 2;
  _thermal_strain_1[_qp] =
      _thermal_strain_1_old[_qp] + (_T[_qp] - _T_old[_qp]) * (alpha1 + alpha1_old) / 2;
  _thermal_strain_2[_qp] =
      _thermal_strain_2_old[_qp] + (_T[_qp] - _T_old[_qp]) * (alpha2 + alpha2_old) / 2;

  _eigenstrain[_qp] = _thermal_strain_0[_qp] * RankTwoTensor::selfOuterProduct(_n0[_qp]) +
                      _thermal_strain_1[_qp] * RankTwoTensor::selfOuterProduct(_n1[_qp]) +
                      _thermal_strain_2[_qp] * RankTwoTensor::selfOuterProduct(_n2[_qp]);
}
