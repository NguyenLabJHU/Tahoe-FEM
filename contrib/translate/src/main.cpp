/* $Id: main.cpp,v 1.4 2001-09-10 20:33:12 sawimme Exp $ */

#include "TranslateIOManager.h"
#include "ExtractIOManager.h"

int main (void)
{
  try 
    {
      int selection;
      cout << "\n1. Datafile Translation \n";
      cout << "2. Nodal Data Extraction to XY Data \n";
      cout << "3. Quadrature Data Extraction (not yet, but soon) \n";
      cout << "\n Select type of translation: ";
      cin >> selection;
      
      TranslateIOManager *dataio;
      StringT program, version;
      switch (selection)
	{
	case 1:
	  {
	    cout << "\n\n Program to translate data files.\n\n";
	    program = "Translate";
	    version = "v1.4";
	    dataio = new TranslateIOManager (cout);
	    break;
	  }
	case 2:
	  {
	    cout << "\n\n Program to extract nodal data.\n\n";
	    program = "Extract";
	    version = "v1.0";
	    dataio = new ExtractIOManager (cout);
	    break;
	  }
	default:
	  throw eGeneralFail;
	}
      dataio->Translate (program, version, program);
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

