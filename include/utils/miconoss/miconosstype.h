#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;

typedef std::vector<double> vecD;
typedef std::vector<std::vector<double>> matrixD;
typedef unsigned int uint;

namespace miconossprint {

inline void printVector(const vecD &V, const std::string &Vname) {
  std::cout << Vname << " = [ ";
  for (uint i = 0; i < V.size(); i++)
    std::cout << V[i] << ", ";
  std::cout << "] \n";
}

inline void printMatrix(const matrixD &M, const std::string &Mname) {
  std::cout << Mname << " = \n[ ";
  for (uint i = 0; i < M.size(); i++) {
    for (uint j = 0; j < M[i].size(); j++)
      std::cout << M[i][j] << ", ";
    if (i + 1 == M.size())
      std::cout << "]\n ";
    else
      std::cout << "\n ";
  }
}

} // namespace miconossprint

namespace miconossmath {

extern "C" void dgesv_(int *, int *, double *, int *, int *, double *, int *,
                       int *);

enum normtype { L2, INF };

inline double L2norm(const vecD &V) {
  double norm = 0;
  for (auto v : V)
    norm += v * v;

  return std::sqrt(norm);
}

inline double LInfnorm(const vecD &V) {
  double norm = 0;
  for (auto v : V) {
    double absv = std::abs(v);
    if (absv > norm)
      norm = absv;
  }
  return norm;
}

inline double norm(const vecD &V, const normtype &nt) {
  double norm = 0;
  switch (nt) {
  case L2:
    norm = L2norm(V);
    break;
  case INF:
    norm = LInfnorm(V);
    break;
  }
  return norm;
}

inline vecD solveAxb(const matrixD &A, const vecD &b, const uint &syssize) {
  int neq = syssize;
  int nrhs = 1;
  vecD b_to_use = b;
  vecD jac(syssize * syssize);
  std::vector<int> IPIV(syssize);
  int info;

  // reorder jacobian in column major
  int k = 0;
  for (uint i = 0; i < syssize; i++)
    for (uint l = 0; l < syssize; l++) {
      jac[k] = A[l][i];
      k += 1;
    }

  dgesv_(&neq, &nrhs, &jac[0], &neq, &IPIV[0], &b_to_use[0], &neq, &info);

  return b_to_use;
}

inline matrixD solveAxNb(const matrixD &A,
                         const std::vector<std::vector<double>> &b,
                         const uint &syssize) {
  int neq = syssize;
  int nrhs = b.size();
  vecD b_to_use(neq * nrhs);
  vecD A_to_use(neq * neq);
  std::vector<int> IPIV(neq);
  int info;

  // reorder jacobian in column major
  int k = 0;
  for (uint i = 0; i < syssize; i++)
    for (uint l = 0; l < syssize; l++) {
      A_to_use[k] = A[l][i];
      k += 1;
    }

  // flattened vector
  k = 0;
  for (uint i = 0; i < nrhs; i++)
    for (uint l = 0; l < syssize; l++) {
      b_to_use[k] = b[i][l];
      k += 1;
    }

  dgesv_(&neq, &nrhs, &A_to_use[0], &neq, &IPIV[0], &b_to_use[0], &neq, &info);

  std::vector<std::vector<double>> b_out(nrhs, std::vector<double>(neq));

  // copy flat solution in to Matrix
  k = 0;
  for (uint i = 0; i < nrhs; i++)
    for (uint l = 0; l < syssize; l++) {
      b_out[i][l] = b_to_use[k];
      k += 1;
    }

  return b_out;
}

} // namespace miconossmath
