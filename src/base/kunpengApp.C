#include "kunpengApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
kunpengApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

kunpengApp::kunpengApp(InputParameters parameters) : MooseApp(parameters)
{
  kunpengApp::registerAll(_factory, _action_factory, _syntax);
}

kunpengApp::~kunpengApp() {}

void
kunpengApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"kunpengApp"});
  Registry::registerActionsTo(af, {"kunpengApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
kunpengApp::registerApps()
{
  registerApp(kunpengApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
kunpengApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  kunpengApp::registerAll(f, af, s);
}
extern "C" void
kunpengApp__registerApps()
{
  kunpengApp::registerApps();
}
