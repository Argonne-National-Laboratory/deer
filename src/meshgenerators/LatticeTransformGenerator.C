//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LatticeTransformGenerator.h"

registerMooseObject("DeerApp", LatticeTransformGenerator);

InputParameters
LatticeTransformGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh to modify");
  params.addClassDescription("Transforms a unit cube mesh to the primitive cell of a Bravais lattice");

  params.addRequiredParam<Point>("a1", "First lattice vector");
  params.addParam<Point>("a2", Point(0,1,0), "Second lattice vector");
  params.addParam<Point>("a3", Point(0,0,1), "Third lattice vector");

  return params;
}

LatticeTransformGenerator::LatticeTransformGenerator(const InputParameters & parameters) :
    MeshGenerator(parameters),
    _a1(getParam<Point>("a1")),
    _a2(getParam<Point>("a2")),
    _a3(getParam<Point>("a3")),
    _input(getMesh("input"))
{

}

std::unique_ptr<MeshBase>
LatticeTransformGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  // Setup transforms
  setupTransforms();

  // Apply to each point
  for (auto & node : mesh->node_ptr_range())
  {
    *node = _T * (*node) + _b;
  }

  return mesh;
}

void
LatticeTransformGenerator::setupTransforms()
{
  // Setup transformation from lattice primitive cell to cube
  Eigen::MatrixXf M(4,4);
  Eigen::MatrixXf V {{0,0,0},{0,0,1},{0,1,1},{1,1,1}};
  
  for (unsigned int i = 0; i < 4; i++)
  {
    Point c = V(i,0) * _a1 + V(i,1) * _a2 + V(i,2) * _a3;
    for (unsigned int j = 0; j < 3; j++)
      M(j,i) = c(j);
  }
  for (unsigned int j = 0; j < 4; j++)
    M(3,j) = 1.0;
  
  Eigen::MatrixXf T = V.transpose() * M.inverse();
  
  RankTwoTensor F;
  Point s;
  for (unsigned int i = 0; i < 3; i++)
  {
    s(i) = T(i,3);
    for (unsigned int j = 0; j < 3; j++)
      F(i,j) = T(i, j);
  }

  // The actual transforms we want are the inverses of these:
  // the maps from the cube to the lattice primitive cell
  _T = F.inverse();
  _b = -_T * s;

}
