#ifndef TIME2DINTERP_H
#define TIME2DINTERP_H

#include "Function.h"

#include "nanoflann.hpp"
#include "KDTreeVectorOfVectorsAdaptor.h"

using namespace nanoflann;

typedef KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<Real>>, Real> kd_tree_t;

class TimeNDInterp;

template <>
InputParameters validParams<TimeNDInterp>();

class TimeNDInterp : public Function
{
 public:
  TimeNDInterp(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) override;

 private:
  /// Read in data from HDF5 file
  void read_data();
  /// Setup the KD Tree
  void setup_kdtree();

 protected:
  size_t ntime_, npoints_, ndim_;
  FileName fname_;
  std::vector<std::vector<Real>> points_;
  std::vector<Real> times_;
  std::vector<std::vector<Real>> values_;
  std::unique_ptr<kd_tree_t> kd_;
};


#endif // TIME2DINTERP_H
