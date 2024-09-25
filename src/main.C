#include "DeerApp.h"
#include "MooseMain.h"

// Begin the main program.
int
main(int argc, char * argv[])
{
  Moose::main<DeerApp>(argc, argv);

  return 0;
}
