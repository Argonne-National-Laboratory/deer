#pragma once

#include "AuxKernel.h"

class IPFColoring : public AuxKernel
{
public:
  static InputParameters validParams();

  IPFColoring(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  /// quaternion
  /// @{
  const VariableValue & _q1;
  const VariableValue & _q2;
  const VariableValue & _q3;
  const VariableValue & _q4;
  /// @}

  /// Sample direction
  const Point _sample_direction;

  /// Crystal symmetry code
  const std::string _crystal_symmetry;

  /// Color component (0=R, 1=G, 2=B)
  const unsigned int _comp;

  /// Sample symmetry code, defaults to "222"
  const std::string _sample_symmetry;

  /// Crystal directions defining triangle
  const Point _v1, _v2, _v3;
};

/// Color points based on quaternion components
Point rgb_ipf_quaternion(Real q1,
                         Real q2,
                         Real q3,
                         Real q4,
                         Point sd,
                         std::string crystal_sym,
                         std::string sample_sym,
                         Point v1,
                         Point v2,
                         Point v3);

/// Actually do the IPF projection and provide a list of points
std::vector<Point>
project_ipf(std::vector<Real> q, Point sd, std::string crystal_sym, std::string sample_sym);

/// Filter a list of points to those in the provided stereographic triangle
std::vector<Point> filter_ipf_points_triangle(std::vector<Point> pts, Point v1, Point v2, Point v3);
