/* $Id: VTKConsoleT.cpp,v 1.1 2001-09-19 01:12:35 recampb Exp $ */

#include "VTKConsoleT.h"

VTKConsoleT::VTKConsoleT(void)
{
  /* set console name */
  iSetName("vtk");

  /* add variables to the console */
  iAddVariable("an_integer", fInteger);
  iAddVariable("a_double", fDouble);  
  iAddVariable("a_string", fString);

  /* add console commands */
  iAddCommand("Integer_Print");
  iAddCommand("String_Print");
}

/* execute given command - returns false on fail */
bool VTKConsoleT::iDoCommand(const StringT& command, StringT& line)
{
  if (command == "Integer_Print")
    {
      cout << "int = " << fInteger << endl;
      return true;
    }
  else if (command == "String_Print")
    {
      cout << "string = " << fString << endl;
      return true;
    }
  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}
