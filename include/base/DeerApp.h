#pragma once

#include "MooseApp.h"

class DeerApp : public MooseApp {
public:
  static InputParameters validParams();
  DeerApp(InputParameters parameters);
  virtual ~DeerApp();

  static void registerApps();
  static void registerAll(Factory &f, ActionFactory &af, Syntax &s);
};
