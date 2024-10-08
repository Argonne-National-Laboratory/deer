#include "AppFactory.h"
#include "DeerApp.h"
#include "ModulesApp.h"
#include "Moose.h"
#include "MooseSyntax.h"

InputParameters
DeerApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

DeerApp::DeerApp(InputParameters parameters) : MooseApp(parameters)
{
  DeerApp::registerAll(_factory, _action_factory, _syntax);
}

DeerApp::~DeerApp() {}

static void
associateSyntaxInner(Syntax & syntax, ActionFactory & /*action_factory*/)
{
  registerSyntax("NEMLMechanicsAction", "NEMLMechanics");
  registerSyntax("RankTwoTensorIntegralAction", "RankTwoTensorIntegralOnDomain/*");
  registerSyntax("RankTwoTensorPostprocessorTimeIntegralAction",
                 "RankTwoTensorPostprocessorTimeIntegral/*");
  registerSyntax("RankTwoTensorPostprocessorTimeDerivativeAction",
                 "RankTwoTensorPostprocessorTimeDerivative/*");
  registerSyntax("CZMStrainAction", "CZMStrain");
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

  associateSyntaxInner(s, af);

  ModulesApp::registerAllObjects<DeerApp>(f, af, s);
}

// External entry point for dynamic application loading
extern "C" void
DeerApp__registerApps()
{
  DeerApp::registerApps();
}

extern "C" void
DeepApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  DeerApp::registerAll(f, af, s);
}
