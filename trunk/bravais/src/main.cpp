// DEVELOPMENT
/* $Id: main.cpp,v 1.6 2002-11-14 01:47:32 saubry Exp $ */
#include <iostream>
#include "ExceptionCodes.h"
#include "PeriodicTableT.h"
#include "StringT.h"

#include "MakeCrystalT.h"

#include "MeshAtom.h"


int main()
{
  MakeCrystalT MC;
  try
    {
      MC.Run();
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


  return 0;
}	
