// $Id: main.cpp,v 1.8 2002-10-28 21:36:33 sawimme Exp $
// created: 6 Oct 1999 by S. A. Wimmer
// program reads input file, runs MakeCSE, writes output file

#include "ExceptionT.h"
#include "MakeCSE_ExecutionT.h"

using namespace Tahoe;

int main (void)
{
  try
    {
      MakeCSE_ExecutionT exe;
      exe.Run ();
      cout << "\n Program Complete\n\n" << endl;
    }

  catch (int ErrorCode) 
    {
      switch (ErrorCode)
	{
	case ExceptionT::kBadInputValue:
	  cout << "\n\n Exiting due to bad input value.\n\n";
	  break;
	case ExceptionT::kOutOfRange:
	  cout << "\n\n Exiting due to index out of range.\n\n";
	  break;
	case ExceptionT::kSizeMismatch:
	  cout << "\n\n Exiting due to a size mismatch.\n\n";
	  break;
	case ExceptionT::kOutOfMemory:
	  cout << "\n\n Out of memory.\n\n";
	  break;
	case ExceptionT::kDatabaseFail:
	  cout << "\n\n Database Error.\n\n";
	  break;
	default:
	  cout << "\n\n Exiting due to failure.\n\n";
	}
      cout << "      Game Over\n\n";
    }
  return 1;
}
