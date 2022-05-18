//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceUserObject.h"

/**
 *  An interface userobject reading boundary based materail property from a
 * text file and storing them for later use
 */
class GBCavitationBoundaryPropertyUO : public InterfaceUserObject
{
public:
  static InputParameters validParams();
  GBCavitationBoundaryPropertyUO(const InputParameters & parameters);
  virtual void initialSetup() override;
  virtual void initialize() override{};
  virtual void execute() override{};
  virtual void finalize() override{};
  virtual void threadJoin(const UserObject & /*uo*/) override{};

  std::map<std::string, Real> getPropertyMap(const dof_id_type elem_id,
                                             const unsigned int side) const;

protected:
  FileName _file_name;
  std::map<std::pair<SubdomainID, SubdomainID>, std::map<std::string, Real>> gb_property_map;
  std::map<std::pair<dof_id_type, unsigned int>, std::pair<SubdomainID, SubdomainID>>
      _elmeside_gbpair_map;
  void ReadGBPropertyFile();
  void constructGBPairsMap();
};
