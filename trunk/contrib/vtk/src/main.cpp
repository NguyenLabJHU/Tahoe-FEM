/* $Id: main.cpp,v 1.1 2001-09-19 01:12:35 recampb Exp $ */

#include "VTKConsoleT.h"
#include "iConsoleT.h"

int main (void)
{
  /* construct VTK console object */
  VTKConsoleT vtk_console, another_console;
  another_console.iSetName("inner_console");
  vtk_console.iAddSub(another_console);

  /* open console */
  iConsoleT("vtk_console.log", vtk_console);
}

