#include "IPFColoring.h"

#include "cp/crystallography.h"
#include "math/rotations.h"
#include "math/tensors.h"

registerMooseObject("DeerApp", IPFColoring);

InputParameters
IPFColoring::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addCoupledVar("q1", "First quaternion component");
  params.addCoupledVar("q2", "Second quaternion component");
  params.addCoupledVar("q3", "Third quaternion component");
  params.addCoupledVar("q4", "Fourth quaternion component");

  params.addRequiredParam<Point>("sample_direction", "Sample direction to project");
  params.addRequiredParam<std::string>("crystal_symmetry", "Crystal symmetry Orbifold code");

  params.addRequiredParam<unsigned int>("component", "Which component of RGB color");

  params.addParam<std::string>("sample_symmetry", "222", "Sample symmetry Orbifold code");

  params.addParam<Point>("v1", Point(0,0,1), "First point of triangle");
  params.addParam<Point>("v2", Point(1,0,1), "Second point of triangle");
  params.addParam<Point>("v3", Point(1,1,1), "Third point of triangle");

  return params;
}

IPFColoring::IPFColoring(const InputParameters & parameters)
  : AuxKernel(parameters),
    _q1(coupledValue("q1")),
    _q2(coupledValue("q2")),
    _q3(coupledValue("q3")),
    _q4(coupledValue("q4")),
    _sample_direction(getParam<Point>("sample_direction")),
    _crystal_symmetry(getParam<std::string>("crystal_symmetry")),
    _comp(getParam<unsigned int>("component")),    
    _sample_symmetry(getParam<std::string>("sample_symmetry")),
    _v1(getParam<Point>("v1")),
    _v2(getParam<Point>("v2")),
    _v3(getParam<Point>("v3"))
{


}

Real
IPFColoring::computeValue()
{
  Point rgb = rgb_ipf_quaternion(_q1[0], _q2[0], _q3[0], _q4[0], _sample_direction,
                                 _crystal_symmetry, _sample_symmetry, _v1, _v2, _v3);

  return rgb(_comp);
}

Point
rgb_ipf_quaternion(Real q1, Real q2, Real q3, Real q4, Point sd, std::string crystal_sym,
                   std::string sample_sym, Point v1, Point v2, Point v3)
{
  auto pts = project_ipf({q1, q2, q3, q4}, sd, crystal_sym, sample_sym);
  auto fpts = filter_ipf_points_triangle(pts, v1, v2, v3);

  // This really shouldn't fail at least for nice crystals...
  if (fpts.size() == 0)
    mooseError("Failed to project an inverse pole figure point");
  
  Point pt = fpts[0];

  // Actually do the color projection

  // Normalize just in case
  v1 /= v1.norm();
  v2 /= v2.norm();
  v3 /= v3.norm();

  std::vector<Point> vs = {v1, v2, v3};
  Point color;
  for (size_t i = 0; i < 3; i++)
    color(i) = vs[i].contract(pt);

  for (size_t i = 0; i < 3; i++) {
    Real mf = 1.0;
    for (size_t j = 0; j < 3; j++) {
      Real t = vs[i].contract(vs[j]);
      if (t < mf)
        mf = t;
    }
    color(i) -= mf;
    color(i) /= (1-mf);
  }

  return color;
}

std::vector<Point>
project_ipf(std::vector<Real> q, Point sd, std::string crystal_sym, std::string sample_sym)
{
  // Setup sample direction
  neml::Vector d({sd(0), sd(1), sd(2)});
  d.normalize();
  
  // Setup orientation
  neml::Orientation qa(q);
  neml::Orientation qo = qa.inverse(); // The following assumes a passive rotation

  // Get crystal and sample symmetry operations
  auto crystal_ops = neml::symmetry_rotations(crystal_sym);
  auto sample_ops = neml::symmetry_rotations(sample_sym);

  std::vector<Point> results;

  for (auto srot : sample_ops) {
    auto cpt = qo.apply(srot.apply(d));
    for (auto crot : crystal_ops) {
      auto fpt = crot.apply(cpt);
      if (fpt(2) > 0)
        results.push_back(Point(fpt(0), fpt(1), fpt(2)));
    }
  }

  return results;
}

std::vector<Point>
filter_ipf_points_triangle(std::vector<Point> pts, Point v1, Point v2, Point v3)
{
  // Normalize just in case
  v1 /= v1.norm();
  v2 /= v2.norm();
  v3 /= v3.norm();
  
  // Get the adjoints
  auto n1 = v1.cross(v2);
  auto n2 = v2.cross(v3);
  auto n3 = v3.cross(v1);

  // Filter
  std::vector<Point> res;
  for (auto pt : pts) {
    if ((pt.contract(n1) >= 0) && (pt.contract(n2) >= 0) && (pt.contract(n3) >= 0))
      res.push_back(pt);
  }

  return res;
}
