#include "Newton.h"

Newton::Newton(NLSystem * _nlsys,
               NLSystemVars * sys_vars,
               const double tolerance,
               const uint _max_iter,
               const miconossmath::normtype normtype)
  : _nlsys(_nlsys),
    _sys_vars(sys_vars),
    _tolerance(tolerance),
    _max_iter(_max_iter),
    _n_eq(_nlsys->getNx()),
    _n_lm(_nlsys->getNmultipliers()),
    _n_total(_n_eq + _n_lm),
    _dx(_n_total),
    _normtype(normtype)
{
}

int
Newton::computeNewtonStep(const vecD & R, const matrixD & J)
{
  int ierr = miconossmath::solveAxb(J, R, _n_total, _dx);

  for (uint i = 0; i < _n_total; i++)
    _dx[i] *= -1;
  return ierr;
}

int
Newton::customSubstepInterruption(NLSystemParameters * const sysparams,
                                  bool & custom_interruption_flag)
{
  custom_interruption_flag = false;
  return 0;
}

int
Newton::solve(vecD & lm, matrixD & J, bool & converged, const bool auto_scale_equation)
{
  converged = false;
  double err = 1e6;
  int ierr = 0;
  uint it = 0;
  if (auto_scale_equation)
    _nlsys->updateEquationScaling();

  _nlsys->updateEquationConstants();
  vecD R = _nlsys->assembleR(lm);
  ierr = miconossmath::norm(R, _normtype, err);
  if (ierr == 0)
  {
    // miconossprint::printVector(R, "intiial residual");
    while (err > _tolerance && it < _max_iter)
    {
      J = _nlsys->assembleJ(lm);
      // miconossprint::printMatrix(J, "Jacobian");
      ierr = computeNewtonStep(R, J);
      if (ierr != 0)
      {
        // std::cerr << "error solving substep  \n";
        break;
      }
      // miconossprint::printVector(_dx, "dx");
      for (uint i = 0; i < _n_eq; i++)
        _sys_vars->setValueFromScaled(i, _sys_vars->getValueScaled(i) + _dx[i]);

      for (uint i = 0; i < _n_lm; i++)
        lm[i] += _dx[_n_eq + i];

      R = _nlsys->assembleR(lm);
      ierr = miconossmath::norm(R, _normtype, err);
      if (ierr != 0)
      {
        // std::cerr << "error is non finite: " << err << " \n";
        break;
      }
      it++;
    }
    if (err < _tolerance && ierr == 0)
    {
      converged = true;
      J = _nlsys->assembleJ(lm);
    }
  }
  return ierr;
}

int
Newton::solveSubstep(vecD & lm,
                     matrixD & J,
                     bool & converged,
                     NLSystemParameters * const sysparams,
                     const std::vector<std::string> & pname,
                     matrixD & Tangent,
                     bool & custom_interruption,
                     double & increment_at_custom_interruption,
                     const uint max_ncut,
                     const bool auto_scale_equation,
                     const bool force_substep)
{
  int ierr = 0;
  converged = false;

  const vecD initVarValues = _sys_vars->getValueVectorOld();
  const double total_increment = sysparams->getValue("dt");
  uint n_cut = 0;
  uint n_total_step = 1;
  if (force_substep)
  {
    n_cut += 1;
    n_total_step *= 2;
    sysparams->setValue("dt", sysparams->getValue("dt") / 2.);
  }
  const uint np_tangent = pname.size();
  do
  {

    const double alpha = 1. / n_total_step;
    _sys_vars->setOldFromVector(initVarValues);
    _sys_vars->setToOld();
    for (uint i = 0; i < _n_lm; i++)
      lm[i] = 0;

    matrixD Tangent_old(np_tangent, vecD(_n_total, 0));
    Tangent = Tangent_old;

    increment_at_custom_interruption = 0;
    ierr = 0;
    for (uint s = 0; s < n_total_step; s++)
    {
      sysparams->setValue("dt_accum", sysparams->getValue("dt") * (s + 1));
      const bool last_substep = s == (n_total_step - 1);
      bool substep_converged = false;
      ierr = solve(lm, J, substep_converged, auto_scale_equation);

      if (substep_converged && ierr == 0)
      {
        if (!last_substep)
          increment_at_custom_interruption += sysparams->getValue("dt");
        else
          increment_at_custom_interruption = total_increment;

        // update global tangent
        Tangent_old = Tangent;
        J = _nlsys->unscaleJacobian(J);
        matrixD dRdP = _nlsys->getDResidualDParams(pname);

        for (uint p = 0; p < np_tangent; p++)
        {
          if (sysparams->isRateParam(pname[p]))
          {
            double dP_dinc = sysparams->getDRateDIncrement(pname[p]);
            for (uint j = 0; j < _n_total; j++)
            {
              dRdP[p][j] *= dP_dinc;
              if (!std::isfinite(dRdP[p][j]))
                std::cerr << "dRdP[p][j] is not finite" << std::to_string(dRdP[p][j]) << " \n ";
            }
          }
        }

        ierr = miconossmath::updateConsistenTangent(J, Tangent_old, dRdP, _n_total, Tangent, alpha);
        if (ierr != 0)
        {
          converged = false;
          std::cerr << "error computing consitent tangent  \n";
          break;
        }

        ierr = customSubstepInterruption(sysparams, custom_interruption);
        if (ierr != 0)
        {
          converged = false;
          std::cerr << "error custom interruption  \n";
          break;
        }

        // if we are at the last step or we receive a custom interruption
        if (last_substep || custom_interruption)
        {
          converged = true;
          break;
        }
        else
        {
          // advance to next step
          _sys_vars->updateOldToCurrent();
          for (uint i = 0; i < _n_lm; i++)
            lm[i] = 0;
        }
      }
      else
      {
        // not converged, break the loop and restart
        converged = false;
        break;
      }
    }
    if (!converged)
    {
      n_cut += 1;
      n_total_step *= 2;
      sysparams->setValue("dt", sysparams->getValue("dt") / 2.);
    }
  } while (n_cut < max_ncut && !converged);

  if (!converged)
  {
    ierr = 1;
    std::cerr << "did not converged after " << std::to_string(n_cut) << " cutabck" << std::endl;
  }
  return ierr;
}
