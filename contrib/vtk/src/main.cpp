/* $Id: main.cpp,v 1.2 2001-09-20 23:11:42 recampb Exp $ */

#include "VTKConsoleT.h"
#include "iConsoleT.h"

int main (void)
{
  /* construct VTK console object */
  VTKConsoleT vtk_console, another_console, color_map;
  another_console.iSetName("inner_console");
  vtk_console.iAddSub(another_console);
  color_map.iSetName("color_map");
  vtk_console.iAddSub(color_map);


  /* open console */
  iConsoleT("vtk_console.log", vtk_console);
}

