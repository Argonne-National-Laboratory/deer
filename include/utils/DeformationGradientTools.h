#pragma once

#include "RankFourTensor.h"
#include "RankTwoTensor.h"

namespace DeformationGradientTools {

RankFourTensor computedRdF(const RankTwoTensor &R, const RankTwoTensor &U) {
  const RankTwoTensor Uhat = U.trace() * RankTwoTensor::Identity() - U;
  unsigned int k, l, m, n, p, q;
  const Real Uhat_det = Uhat.det();

  RankFourTensor dR_dF;
  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++)
      for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++) {
          dR_dF(k, l, m, n) = 0.;
          for (p = 0; p < 3; p++)
            for (q = 0; q < 3; q++)
              dR_dF(k, l, m, n) +=
                  R(k, p) * (Uhat(p, q) * R(m, q) * Uhat(n, l) -
                             Uhat(p, n) * R(m, q) * Uhat(q, l));

          dR_dF(k, l, m, n) /= Uhat_det;
        }

  return dR_dF;
}
RankFourTensor computedFinversedF(const RankTwoTensor &F_inv) {
  return -F_inv.mixedProductIkJl(F_inv.transpose());
};

} // namespace DeformationGradientTools
