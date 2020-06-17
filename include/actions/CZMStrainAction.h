#pragma once

#include "Action.h"

class CZMStrainAction : public Action {
public:
  static InputParameters validParams();
  CZMStrainAction(const InputParameters &params);

  virtual void act();

protected:
  void addInterfaceStrainMaterial();
  void computeScalingVolume();
  void addInterfaceStrainRateAction();
  void addIntegrateInterfaceStrainRateAction();
  void addEquivalentStrain(const PostprocessorName &rank_two_base_name);

  const std::vector<VariableName> _displacements;
  const std::vector<SubdomainName> _block;
  const std::vector<BoundaryName> _boundary;
  const bool _scaled;
  const PostprocessorName _bulk_volume_PP;
  const PostprocessorName _czm_strain_base_name;
  const bool _compute_cumulative_strain;
  const bool _compute_equivalent_strain;
  const bool _ld;
};
