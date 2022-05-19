#pragma once

#include "Action.h"

class CZMStrainAction : public Action
{
public:
  static InputParameters validParams();
  CZMStrainAction(const InputParameters & params);

  virtual void act();

protected:
  void addInterfaceStrainMaterial();
  void computeScalingVolume();
  void addInterfaceStrain();
  void addInterfaceStrainRate();
  void addEquivalentStrain(const PostprocessorName & rank_two_base_name);

  const std::vector<SubdomainName> _block;
  const std::vector<BoundaryName> _boundary;
  const bool _scaled;
  const PostprocessorName _bulk_volume_PP;
  const PostprocessorName _area_ratio_PP;
  const PostprocessorName _czm_strain_scale_PP;
  const PostprocessorName _czm_strain_base_name;
  const bool _compute_czm_strain_rate;
  const bool _compute_equivalent_strain;

  // map between tensor components and names
  const std::map<std::pair<int, int>, std::string> _tensor_map = {{std::make_pair(0, 0), "xx"},
                                                                  {std::make_pair(1, 1), "yy"},
                                                                  {std::make_pair(2, 2), "zz"},
                                                                  {std::make_pair(0, 1), "xy"},
                                                                  {std::make_pair(0, 2), "xz"},
                                                                  {std::make_pair(1, 2), "yz"}};

  std::vector<MaterialPropertyName> _czm_mp_strain_names = {
      "czm_total_strain", "czm_normal_strain", "czm_sliding_strain"};
};
