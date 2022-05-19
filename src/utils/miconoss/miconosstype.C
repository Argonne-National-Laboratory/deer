#pragma once

#include "miconosstype.h"

namespace miconossprint
{

void
printVector(const vecD & V, const std::string & Vname)
{
  std::cerr << Vname << " = [ ";
  for (uint i = 0; i < V.size(); i++)
    std::cerr << V[i] << ", ";
  std::cerr << "] \n";
}

void
printMatrix(const matrixD & M, const std::string & Mname)
{
  std::cerr << Mname << " = \n[ ";
  for (uint i = 0; i < M.size(); i++)
  {
    for (uint j = 0; j < M[i].size(); j++)
      std::cerr << M[i][j] << ", ";
    if (i + 1 == M.size())
      std::cerr << "]\n ";
    else
      std::cerr << "\n ";
  }
}

} // namespace miconossprint

namespace miconossmath
{

int
L2norm(const vecD & V, double & norm)
{
  int ierr = 0;
  norm = 0;
  for (auto v : V)
  {
    if (std::isfinite(v))
      norm += v * v;
    else
    {
      ierr = 1;
      norm = std::numeric_limits<double>::quiet_NaN();
      break;
    }
  }
  if (ierr == 0)
    norm = std::sqrt(norm);

  return ierr;
}

int
LInfnorm(const vecD & V, double & norm)
{
  norm = 0;
  int ierr = 0;
  for (auto v : V)
  {
    if (std::isfinite(v))
    {
      double absv = std::abs(v);
      if (absv > norm)
        norm = absv;
    }
    else
    {
      norm = std::numeric_limits<double>::quiet_NaN();
      ierr = 1;
      break;
    }
  }
  return ierr;
}

int
norm(const vecD & V, const normtype & nt, double & norm)
{
  norm = 0;
  int ierr = 0;
  switch (nt)
  {
    case L2:
      ierr = L2norm(V, norm);
      break;
    case INF:
      ierr = LInfnorm(V, norm);
      break;
  }
  return ierr;
}

int
solveAxb(const matrixD & A, const vecD & b, const uint & syssize, vecD & b_out)
{
  int neq = syssize;
  int nrhs = 1;
  b_out = b;
  vecD jac(syssize * syssize);
  std::vector<int> IPIV(syssize);
  int info;

  // reorder jacobian in column major
  int k = 0;
  for (uint i = 0; i < syssize; i++)
    for (uint l = 0; l < syssize; l++)
    {
      jac[k] = A[l][i];
      k += 1;
    }

  dgesv_(&neq, &nrhs, &jac[0], &neq, &IPIV[0], &b_out[0], &neq, &info);

  if (info < 0)
    std::cerr << "solveAxb the " + std::to_string(info) + "the argument has an illegal value";

  if (info > 0)
    std::cerr << "solveAxb factorization U si singular \n";

  return info;
}

int
solveAxNb(const matrixD & A, const matrixD & b, const uint & syssize, matrixD & b_out)
{
  int neq = syssize;
  int nrhs = b.size();
  vecD b_to_use(neq * nrhs);
  vecD A_to_use(neq * neq);
  std::vector<int> IPIV(neq);
  int info;

  // reorder jacobian in column major
  int k = 0;
  for (uint i = 0; i < syssize; i++)
    for (uint l = 0; l < syssize; l++)
    {
      A_to_use[k] = A[l][i];
      k += 1;
    }

  // flattened vector
  k = 0;
  for (int i = 0; i < nrhs; i++)
    for (uint l = 0; l < syssize; l++)
    {
      b_to_use[k] = b[i][l];
      k += 1;
    }

  dgesv_(&neq, &nrhs, &A_to_use[0], &neq, &IPIV[0], &b_to_use[0], &neq, &info);

  if (info < 0)
    std::cerr << "solveAxNb the " + std::to_string(-info) + "th argument has an illegal value";

  if (info > 0)
    std::cerr << "solveAxNb factorization U si singular \n";

  if (info == 0)
  {
    b_out = matrixD(nrhs, vecD(neq));

    // copy flat solution in to Matrix
    k = 0;
    for (int i = 0; i < nrhs; i++)
      for (uint l = 0; l < syssize; l++)
      {
        b_out[i][l] = b_to_use[k];
        k += 1;
      }
  }
  return info;
}

int
updateConsistenTangent(const matrixD & J,
                       matrixD & TangentOld,
                       matrixD & dRdP,
                       const uint & syssize,
                       matrixD & NewTangent,
                       const double alpha)
{

  const uint nparam = dRdP.size();

  for (uint p = 0; p < nparam; p++)
    for (uint j = 0; j < syssize; j++)
      TangentOld[p][j] -= alpha * dRdP[p][j];

  int ierr = solveAxNb(J, TangentOld, syssize, NewTangent);

  for (uint p = 0; p < nparam; p++)
    for (uint j = 0; j < syssize; j++)
      if (!std::isfinite(NewTangent[p][j]))
      {
        ierr = 1;
        break;
      }

  return ierr;
}

} // namespace miconossmath
