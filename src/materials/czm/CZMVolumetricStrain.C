#include "CZMVolumetricStrain.h"
#include "CohesiveZoneModelTools.h"
#include "RotationMatrix.h"

registerMooseObject("DeerApp", CZMVolumetricStrain);

InputParameters CZMVolumetricStrain::validParams() {
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription("Material computing volumetric strain resulting "
                             "from interface separation and sliding."
                             "It also provides the normal and sliding strain "
                             "contrbuting to the total interface strain.");
  MooseEnum strainType("SMALL FINITE", "SMALL");
  params.addParam<MooseEnum>("strain", strainType, "Strain formulation");
  params.addParam<std::string>("base_name", "Material property base name");
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

CZMVolumetricStrain::CZMVolumetricStrain(const InputParameters &parameters)
    : InterfaceMaterial(parameters), _normals(_assembly.normals()),
      _base_name(isParamValid("base_name") &&
                         !getParam<std::string>("base_name").empty()
                     ? getParam<std::string>("base_name") + "_"
                     : ""),
      _displacement_jump_global(getMaterialPropertyByName<RealVectorValue>(
          _base_name + "displacement_jump_global")),
      _czm_total_rotation(getMaterialPropertyByName<RankTwoTensor>(
          _base_name + "czm_total_rotation")),
      _czm_total_strain(declarePropertyByName<RankTwoTensor>(
          _base_name + "czm_total_strain")),
      _czm_normal_strain(declarePropertyByName<RankTwoTensor>(
          _base_name + "czm_normal_strain")),
      _czm_sliding_strain(declarePropertyByName<RankTwoTensor>(
          _base_name + "czm_sliding_strain")),
      _dadA_mp(declarePropertyByName<Real>(_base_name + "czm_area_ratio")),
      _strain(getParam<MooseEnum>("strain").getEnum<Strain>()),
      _F_czm(
          _strain == Strain::Finite
              ? &getMaterialPropertyByName<RankTwoTensor>(_base_name + "F_czm")
              : nullptr) {}

void CZMVolumetricStrain::initQpStatefulProperties() {
  _czm_total_strain[_qp] = 0;
  _czm_normal_strain[_qp] = 0;
  _czm_sliding_strain[_qp] = 0;
}

void CZMVolumetricStrain::computeQpProperties() {
  // initialize kinematics variable for small deformation
  _dadA_mp[_qp] = 1;

  if (_strain == Strain::Finite)
    _dadA_mp[_qp] = CohesiveZoneModelTools::computeAreaRatio(
        (*_F_czm)[_qp].inverse().transpose(), (*_F_czm)[_qp].det(),
        _normals[_qp]);

  computeInterfaceStrain();
}

void CZMVolumetricStrain::computeInterfaceStrain() {
  const RealVectorValue n_average =
      _czm_total_rotation[_qp] * RealVectorValue(1.0, 0.0, 0.0);
  const RankTwoTensor u_outer_n =
      outer_product(_displacement_jump_global[_qp], n_average);
  const RankTwoTensor n_outer_n = outer_product(n_average, n_average);

  _czm_total_strain[_qp] =
      (u_outer_n + u_outer_n.transpose()) / 2. * _dadA_mp[_qp];
  _czm_normal_strain[_qp] = _displacement_jump_global[_qp] * n_average *
                            (n_outer_n + n_outer_n.transpose()) / 2. *
                            _dadA_mp[_qp];
  _czm_sliding_strain[_qp] = _czm_total_strain[_qp] - _czm_normal_strain[_qp];
}
