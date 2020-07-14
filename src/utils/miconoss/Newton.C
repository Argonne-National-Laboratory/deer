#include "Newton.h"

Newton::Newton(NLSystem *_nlsys, NLSystemVars *sys_vars, const double tolerance,
               const uint _max_iter, const miconossmath::normtype normtype)
    : _nlsys(_nlsys), _sys_vars(sys_vars), _tolerance(tolerance),
      _max_iter(_max_iter), _n_eq(_nlsys->getNx()),
      _n_lm(_nlsys->getNmultipliers()), _n_total(_n_eq + _n_lm), _dx(_n_total),
      _normtype(normtype) {}

int Newton::computeNewtonStep(const vecD &R, const matrixD &J) {
  int ierr = miconossmath::solveAxb(J, R, _n_total, _dx);

  for (uint i = 0; i < _n_total; i++)
    _dx[i] = -_dx[i];
  return ierr;
}

bool Newton::solve(vecD &lm, matrixD &J, const bool auto_scale_equation) {
  bool converged = false;
  uint it = 0;
  int ierr = 0;
  if (auto_scale_equation)
    _nlsys->updateEquationScaling();

  _nlsys->updateEquationConstants();
  vecD R = _nlsys->assembleR(lm);
  double err = miconossmath::norm(R, _normtype);
  if (std::isfinite(err)) {
    // miconossprint::printVector(R, "intiial residual");
    while (err > _tolerance && it < _max_iter) {
      J = _nlsys->assembleJ(lm);
      // miconossprint::printMatrix(J, "Jacobian");
      ierr = computeNewtonStep(R, J);
      if (ierr != 0) {
        std::cerr << "error solving substep  \n";
        break;
      }
      // miconossprint::printVector(_dx, "dx");
      for (uint i = 0; i < _n_eq; i++)
        _sys_vars->setValueFromScaled(i, _sys_vars->getValueScaled(i) + _dx[i]);

      for (uint i = 0; i < _n_lm; i++)
        lm[i] += _dx[_n_eq + i];

      R = _nlsys->assembleR(lm);
      err = miconossmath::norm(R, _normtype);
      if (!std::isfinite(err))
        break;
      it++;
    }
    if (std::isfinite(err) && err < _tolerance && ierr == 0) {
      converged = true;
      J = _nlsys->assembleJ(lm);
    }
  }
  return converged;
}

bool Newton::solveSubstep(vecD &lm, matrixD &J,
                          NLSystemParameters *const sysparams,
                          const std::vector<std::string> &pname,
                          matrixD &Tangent, bool &custom_interruption,
                          double &increment_at_custom_interruption,
                          const uint max_ncut, const bool auto_scale_equation,
                          const bool force_substep) {
  bool converged = false;
  const vecD initVarValues = _sys_vars->getValueVectorOld();
  const double total_increment = sysparams->getValue("dt");
  uint n_cut = 0;
  uint n_total_step = 1;
  if (force_substep) {
    n_cut += 1;
    n_total_step *= 2;
    sysparams->setValue("dt", sysparams->getValue("dt") / 2.);
  }
  const uint np_tangent = pname.size();
  do {

    const double alpha = 1. / n_total_step;
    _sys_vars->setOldFromVector(initVarValues);
    _sys_vars->setToOld();
    for (uint i = 0; i < _n_lm; i++)
      lm[i] = 0;

    matrixD Tangent_old(np_tangent, vecD(_n_total, 0));
    Tangent = Tangent_old;

    increment_at_custom_interruption = 0;
    for (uint s = 0; s < n_total_step; s++) {
      const bool last_substep = s == (n_total_step - 1);
      bool substep_converged = solve(lm, J, auto_scale_equation);
      if (!last_substep)
        increment_at_custom_interruption += sysparams->getValue("dt");
      else
        increment_at_custom_interruption = total_increment;

      if (substep_converged) {

        // update global tangent
        Tangent_old = Tangent;
        J = _nlsys->unscaleJacobian(J);
        matrixD dRdP = _nlsys->getDResidualDParams(pname);

        // if we have rate paramters, and we want the consistent tangent be
        // we have to scale back dRdP using the derivative of the rate paramter
        // w.r.t. the increment

        for (uint p = 0; p < np_tangent; p++) {
          double dP_dinc = 1.;
          if (sysparams->isRateParam(pname[p]))
            dP_dinc = sysparams->getDRateDIncrement(pname[p]);
          for (uint j = 0; j < _n_total; j++)
            Tangent_old[p][j] -= alpha * dRdP[p][j] * dP_dinc;
        }

        int ierr = miconossmath::solveAxNb(J, Tangent_old, _n_total, Tangent);
        if (ierr != 0) {
          std::cerr << "error computing consitent tangent  \n";
          break;
        }
        custom_interruption = customSubstepInterruption(sysparams);

        // if we are at the last step or we receive a custom interruption
        if (last_substep || custom_interruption)
          converged = true;
        else {
          // advance to next step
          _sys_vars->updateOldToCurrent();
          for (uint i = 0; i < _n_lm; i++)
            lm[i] = 0;
        }

      } else {
        // not converged, break the loop and restart
        converged = false;
        break;
      }
    }
    if (!converged) {
      n_cut += 1;
      n_total_step *= 2;
      sysparams->setValue("dt", sysparams->getValue("dt") / 2.);
    }
  } while (n_cut < max_ncut && !converged);

  if (!converged)
    std::cerr << "did not converged after " << std::to_string(n_cut)
              << " cutabck" << std::endl;

  return converged;
}
