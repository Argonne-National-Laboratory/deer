#ifndef COMPUTERADIATIONSWELLINGEIGENSTRAIN_H
#define COMPUTERADIATIONSWELLINGEIGENSTRAIN_H

#include "ComputeEigenstrainBase.h"

class ComputeRadiationSwellingEigenstrain;

template <>
InputParameters validParams<ComputeRadiationSwellingEigenstrain>();

/**
 *  Calculate an eigenstrain based on two functions:
 *  1) Dose rate versus time
 *  2) Swelling versus dose
 *
 *  The swelling curve is input as a swelling relation, i.e. 
 *  change in volume / volume and the class takes care of the conversion
 *  to strain
 */
class ComputeRadiationSwellingEigenstrain: public ComputeEigenstrainBase
{
 public:
  ComputeRadiationSwellingEigenstrain(const InputParameters & parameters);
  virtual void initQpStatefulProperties() override;
  virtual void computeQpEigenstrain() override;

 private:
  Function & _swelling;
  Function & _dose_rate;
  MaterialProperty<Real> & _dose;
  const MaterialProperty<Real> & _dose_old;
};


#endif
