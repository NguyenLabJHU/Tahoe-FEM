/* $Id: main.cpp,v 1.5 2001-11-07 02:34:45 paklein Exp $ */

#include "VTKConsoleT.h"
#include "iConsoleT.h"

int main (int argc, char* argv[])
{
  /* list of command-line arguments */
  ArrayT<StringT> arguments(argc);
  for (int i = 0; i < arguments.Length(); i++)
	arguments[i] = argv[i];

  /* construct VTK console object */
  VTKConsoleT vtk_console(arguments);

  /* open console */
  iConsoleT("vtk_console.log", vtk_console);
}

