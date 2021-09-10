#!/bin/bash

# an example of how to generate a 10 grain mesh in neper with equiaxed grains

~/path/to/neper/neper -T -n 10 -id 25 -ori random -oricrysym cubic -oridescriptor euler-kocks -regularization 1 -morpho 'diameq:lognormal(1,0.005),1-sphericity:lognormal(0.145,0.03)'  -o test_10_gg  -domain "cube(0.093,0.093,0.093)"

~/path/to/neper/neper -M -elttype tet -order 1 -format vtk test_10_gg.tess  -o test_10_gg_rcl1pt1 -meshqualmin 1 -rcl 1.1

# after you get the vtk from neper you need to convert it to exodus using neper2warp
