/* $Id: main.cpp,v 1.1 2001-09-01 00:05:37 paklein Exp $ */

#include "ios_fwd_decl.h"
#include "StringT.h"
#include "TranslateIOManager.h"
#include "ifstreamT.h"

void main (void)
{
  try
    {
      // determine input file 
      char *program_name = "tetexo";
      char *program_version = "v1.3 (10 May 2000)";
      TranslateIOManager iodata (cout);
      ifstreamT in ('#');
      bool interactive = true;
      iodata.ReadParameters (in, true, program_name, program_version);
      iodata.Translate ();
      cout << "\nProgram Complete\n\n" << endl;
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
	  //case eGeneralFail:
	  //case eStop:
	  //default:
	}
      cout << "      Game Over\n\n";
    }
}

