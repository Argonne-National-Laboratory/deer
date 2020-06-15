#pragma once

#include "Action.h"

class CZMStrainAction : public Action {
public:
  static InputParameters validParams();
  CZMStrainAction(const InputParameters &params);

  virtual void act();

protected:
  void addInterfaceStrainRateAction();
  void addIntegrateInterfaceStrainRateAction();

  const std::vector<SubdomainName> _block;
  const std::vector<BoundaryName> _boundary;
  const bool _scaled;
  const PostprocessorName _bulk_volume_PP;
  const PostprocessorName _czm_strain_base_name;
  const bool _add_integrated_interface_strains;
};
