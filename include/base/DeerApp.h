#ifndef DEERAPP_H
#define DEERAPP_H

#include "MooseApp.h"

class DeerApp;

template<>
InputParameters validParams<DeerApp>();

class DeerApp : public MooseApp
{
public:
  DeerApp(InputParameters parameters);
  virtual ~DeerApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* DEERAPP_H */
