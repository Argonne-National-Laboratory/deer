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

inline void printVector(const vecD &V, const std::string &Vname) {
  std::cerr << Vname << " = [ ";
  for (uint i = 0; i < V.size(); i++)
    std::cerr << V[i] << ", ";
  std::cerr << "] \n";
}

inline void printMatrix(const matrixD &M, const std::string &Mname) {
  std::cerr << Mname << " = \n[ ";
  for (uint i = 0; i < M.size(); i++) {
    for (uint j = 0; j < M[i].size(); j++)
      std::cerr << M[i][j] << ", ";
    if (i + 1 == M.size())
      std::cerr << "]\n ";
    else
      std::cerr << "\n ";
  }
}

} // namespace miconossprint

namespace miconossmath {

extern "C" void dgesv_(int *, int *, double *, int *, int *, double *, int *,
                       int *);

enum normtype { L2, INF };

inline double L2norm(const vecD &V) {
  double norm = 0;
  for (auto v : V) {
    if (std::isfinite(v))
      norm += v * v;
    else {
      norm = std::numeric_limits<double>::quiet_NaN();
      break;
    }
  }
  if (std::isfinite(norm))
    norm = std::sqrt(norm);

  return norm;
}

inline double LInfnorm(const vecD &V) {
  double norm = 0;
  for (auto v : V) {
    if (std::isfinite(v)) {
      double absv = std::abs(v);
      if (absv > norm)
        norm = absv;
    } else {
      norm = std::numeric_limits<double>::quiet_NaN();
      break;
    }
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

inline int solveAxb(const matrixD &A, const vecD &b, const uint &syssize,
                    vecD &b_out) {
  int neq = syssize;
  int nrhs = 1;
  b_out = b;
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

  dgesv_(&neq, &nrhs, &jac[0], &neq, &IPIV[0], &b_out[0], &neq, &info);

  if (info < 0)
    std::cerr << "solveAxb the " + std::to_string(info) +
                     "th argument has an illegal value";

  if (info > 0)
    std::cerr << "solveAxb factorization U si singular \n";

  return info;
}

inline int solveAxNb(const matrixD &A,
                     const std::vector<std::vector<double>> &b,
                     const uint &syssize, matrixD &b_out) {
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
  for (int i = 0; i < nrhs; i++)
    for (uint l = 0; l < syssize; l++) {
      b_to_use[k] = b[i][l];
      k += 1;
    }

  dgesv_(&neq, &nrhs, &A_to_use[0], &neq, &IPIV[0], &b_to_use[0], &neq, &info);

  if (info < 0)
    std::cerr << "solveAxNb the " + std::to_string(-info) +
                     "th argument has an illegal value";

  if (info > 0)
    std::cerr << "solveAxNb factorization U si singular \n";

  if (info == 0) {
    b_out = matrixD(nrhs, vecD(neq));

    // copy flat solution in to Matrix
    k = 0;
    for (int i = 0; i < nrhs; i++)
      for (uint l = 0; l < syssize; l++) {
        b_out[i][l] = b_to_use[k];
        k += 1;
      }
  }
  return info;
}

} // namespace miconossmath
