/* $Id: main.cpp,v 1.3 2001-09-10 16:57:45 sawimme Exp $ */

#include "TranslateIOManager.h"

int main (void)
{
  try 
    {
      cout << "\n\n Program to translate data files.\n\n";
      StringT program = "Translate";
      StringT version = "v1.4";
      TranslateIOManager dataio (cout);
      dataio.Translate (program, version, program);
      cout << "\n\n Progam Complete.\n\n";
    }
  catch (int ErrorCode)
    {
      cout << "\n\n Exiting due to error . . . ";
      switch (ErrorCode)
	{
	case eBadInputValue:
	  cout << " Bad Input Value\n";
	  break;
	case eOutOfRange:
	  cout << " Out of Range\n";
	  break;
	case eSizeMismatch:
	  cout << " Size Mismatch\n";
	  break;
	case eOutOfMemory:
	  cout << " Out of Memory\n";
	  break;
	case eDatabaseFail:
	  cout << " Error with database\n";
	  break;
	}
      cout << "\n\n Game Over\n\n";
    }
  return 1;
}

