#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>
using namespace std;

/// miconoss vector
typedef std::vector<double> vecD;
/// miconoss Matrix
typedef std::vector<std::vector<double>> matrixD;
/// miconoss unsigned int
typedef unsigned int uint;

namespace miconossprint
{
/// method used to print a VecD
void printVector(const vecD & V, const std::string & Vname);
/// method used to print a matrixD
void printMatrix(const matrixD & M, const std::string & Mname);

} // namespace miconossprint

namespace miconossmath
{

extern "C" void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);

/// the allowed norm types
enum normtype
{
  L2,
  INF
};

/// method computing the L2 norm
int L2norm(const vecD & V, double & norm);
/// method computing the infinty norm
int LInfnorm(const vecD & V, double & norm);

/// compute the norm of avector based on the norm type
int norm(const vecD & V, const normtype & nt, double & norm);

/// method soving Ax=b with be being being a vector
int solveAxb(const matrixD & A, const vecD & b, const uint & syssize, vecD & b_out);

/// method soving Ax=b with be being a non square matrix
int solveAxNb(const matrixD & A, const matrixD & b, const uint & syssize, matrixD & b_out);

/// method used to update the consistent tangent
int updateConsistenTangent(const matrixD & J,
                           matrixD & TangentOld,
                           matrixD & dRdP,
                           const uint & syssize,
                           matrixD & NewTangent,
                           const double alpha);

} // namespace miconossmath
