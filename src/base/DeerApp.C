#include "DeerApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

// New functions
#include "ThicknessGradient.h"
#include "CapGradient.h"

// New materials
#include "ComputeNEMLStress.h"
#include "ComputeThermalExpansionEigenstrainNEML.h"
#include "ComputeRadiationSwellingEigenstrain.h"

template<>
InputParameters validParams<DeerApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

DeerApp::DeerApp(InputParameters parameters) :
    MooseApp(parameters)
{
  DeerApp::registerAll(_factory, _action_factory, _syntax);
}

DeerApp::~DeerApp()
{
}

void
DeerApp::registerApps()
{
  registerApp(DeerApp);
}

void
DeerApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  Registry::registerObjectsTo(f, {"DeerApp"});
  Registry::registerActionsTo(af, {"DeerApp"});

  ModulesApp::registerAll(f, af, s);
}

// External entry point for dynamic application loading
extern "C" void
DeerApp__registerApps()
{
  DeerApp::registerApps();
}

extern "C" void
DeepApp__registerAll(Factory & f, ActionFactory &af, Syntax & s)
{
  DeerApp::registerAll(f, af, s);
}
