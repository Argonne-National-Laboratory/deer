#include "DeerApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

// New functions
#include "TimeNDInterp.h"
#include "ThicknessGradient.h"
#include "CapGradient.h"

// New materials
#include "ComputeNEMLStress.h"
#include "ComputeThermalExpansionEigenstrainNEML.h"

template<>
InputParameters validParams<DeerApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

DeerApp::DeerApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  DeerApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  DeerApp::associateSyntax(_syntax, _action_factory);
}

DeerApp::~DeerApp()
{
}

// External entry point for dynamic application loading
extern "C" void DeerApp__registerApps() { DeerApp::registerApps(); }
void
DeerApp::registerApps()
{
  registerApp(DeerApp);
}

// External entry point for dynamic object registration
extern "C" void DeerApp__registerObjects(Factory & factory) { DeerApp::registerObjects(factory); }
void
DeerApp::registerObjects(Factory & factory)
{
  registerFunction(TimeNDInterp);
  registerFunction(ThicknessGradient);
  registerFunction(CapGradient);
  registerMaterial(ComputeNEMLStress);
  registerMaterial(ComputeThermalExpansionEigenstrainNEML);
}

// External entry point for dynamic syntax association
extern "C" void DeerApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { DeerApp::associateSyntax(syntax, action_factory); }
void
DeerApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
