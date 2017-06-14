#include "TimeNDInterp.h"

#include "hdf5.h"

#include <algorithm>

template <>
InputParameters
validParams<TimeNDInterp>()
{
  InputParameters params = validParams<Function>();
  params.addRequiredParam<FileName>("file", "File containing data.");
  return params;
}

TimeNDInterp::TimeNDInterp(const InputParameters & parameters)
  : Function(parameters), fname_(getParam<FileName>("file"))
{
  read_data();
  setup_kdtree();
}

Real TimeNDInterp::value(Real t, const Point & p)
{
  // Need a vector
  std::vector<Real> query(ndim_);
  // Manual, should fix
  query[0] = p(0);
  query[1] = p(1);

  // Now this is just a call to the kD tree
  std::vector<size_t> indexes(1);
  std::vector<Real> dists(1);

  nanoflann::KNNResultSet<Real> resultSet(1);

  resultSet.init(&indexes[0], &dists[0]);
  kd_->index->findNeighbors(resultSet, &query[0], nanoflann::SearchParams(10));
  
  size_t close = indexes[0];

  // Now figure out the time interpolation
  if (t < times_[0]) {
    return values_[0][close];
  }
  else if (t > times_.back()) {
    return values_.back()[close];
  }
  else {
    auto up = std::upper_bound(times_.begin(), times_.end(), t);
    size_t index = up - times_.begin();
    double t1 = times_[index-1];
    double t2 = times_[index];
    double v1 = values_[index-1][close];
    double v2 = values_[index][close];

    return v1 + (t - t1) * (v2 - v1) / (t2 - t1);
  }
}

void TimeNDInterp::read_data()
{
  // Setup
  herr_t status;
  hid_t file_id;
  
  file_id = H5Fopen(fname_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // Read in the times
  hid_t dset, dspace;
  size_t ndims;
  dset = H5Dopen(file_id, "t", H5P_DEFAULT);
  dspace = H5Dget_space(dset);
  ndims = H5Sget_simple_extent_ndims(dspace);
  assert(ndims == 1);
  hsize_t dims_temp[1];
  H5Sget_simple_extent_dims(dspace, dims_temp, NULL);
  ntime_ = dims_temp[0];
  times_.resize(ntime_);
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
          &(times_[0]));

  // Read in the points
  dset = H5Dopen(file_id, "xy", H5P_DEFAULT);
  dspace = H5Dget_space(dset);
  ndims = H5Sget_simple_extent_ndims(dspace);
  assert(ndims == 2);
  hsize_t dims_pts[2];
  H5Sget_simple_extent_dims(dspace, dims_pts, NULL);
  npoints_ = dims_pts[0];
  ndim_ = dims_pts[1];
  
  // Buffer
  double * buffer;
  buffer = new double[npoints_ * ndim_];

  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  points_.resize(npoints_);
  for (int i=0; i < npoints_; i++) {
    points_[i].resize(ndim_);
    for (int j=0; j < ndim_; j++) {
      points_[i][j] = buffer[i*ndim_+j];
    }
  }

  delete [] buffer;

  // Read in the temperatures
  dset = H5Dopen(file_id, "T", H5P_DEFAULT);
  dspace = H5Dget_space(dset);
  ndims = H5Sget_simple_extent_ndims(dspace);
  assert(ndims == 2);
  hsize_t dims_temps[2];
  H5Sget_simple_extent_dims(dspace, dims_temps, NULL);
  assert(dims_temps[0] == ntime_);
  assert(dims_temps[1] == npoints_);
  
  // Buffer
  buffer = new double[ntime_ * npoints_];
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  values_.resize(ntime_);
  for (int i=0; i< ntime_; i++) {
    values_[i].resize(npoints_);
    for (int j=0; j < npoints_; j++) {
      values_[i][j] = buffer[i*npoints_+j];
    }
  }

  delete [] buffer;

  // Close down
  status = H5Fclose(file_id);
}

void TimeNDInterp::setup_kdtree() {
  // Pretty simple library call
  std::unique_ptr<kd_tree_t> obj(new kd_tree_t(ndim_, points_, 10));
  kd_ = std::move(obj);
  kd_->index->buildIndex();
}
