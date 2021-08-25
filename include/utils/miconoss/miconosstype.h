#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>
using namespace std;

typedef std::vector<double> vecD;
typedef std::vector<std::vector<double>> matrixD;
typedef unsigned int uint;

namespace miconossprint {

void printVector(const vecD &V, const std::string &Vname);

void printMatrix(const matrixD &M, const std::string &Mname);

} // namespace miconossprint

namespace miconossmath {

extern "C" void dgesv_(int *, int *, double *, int *, int *, double *, int *,
                       int *);

enum normtype { L2, INF };

int L2norm(const vecD &V, double &norm);

int LInfnorm(const vecD &V, double &norm);

int norm(const vecD &V, const normtype &nt, double &norm);

int solveAxb(const matrixD &A, const vecD &b, const uint &syssize, vecD &b_out);

int solveAxNb(const matrixD &A, const matrixD &b, const uint &syssize,
              matrixD &b_out);

int updateConsistenTangent(const matrixD &J, matrixD &TangentOld, matrixD &dRdP,
                           const uint &syssize, matrixD &NewTangent,
                           const double alpha);

} // namespace miconossmath
