//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EffectiveStressTools.h"

namespace EffectiveStressTools {

MooseEnum scalarOptions() {
  return MooseEnum("VonMises Hydrostatic Huddleston Hayhurst MaxPrincipal "
                   "Tresca RCCMRXMises RCCMRXTresca maxS1AndMises");
}

} // namespace EffectiveStressTools
