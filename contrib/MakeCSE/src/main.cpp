// $Id: main.cpp,v 1.7 2002-10-25 21:02:59 paklein Exp $
// created: 6 Oct 1999 by S. A. Wimmer
// program reads input file, runs MakeCSE, writes output file

#include "MakeCSE_FEManager.h"
#include "ExceptionCodes.h"
#include "ifstreamT.h"

using namespace Tahoe;

bool IsInteractive (ifstreamT& in);

int main (void)
{
  try
    {
      /* determine input file */
      const StringT program_name = "MakeCSE";
      const StringT program_version = "v5.1 (Oct 2002)";
 
      /* create input/output manager and read input file */
      cout << "\n Welcome to: " << program_name << " " << program_version;
      cout << "\n\n Build Date: " << __DATE__ " " << __TIME__ << "\n\n";

      ifstreamT in ('#');
      bool interactive = IsInteractive (in);

      StringT outfile (81);
      if (interactive)
	{
	  cout << "\n Enter output file: ";
	  cin.getline (outfile.Pointer(), 80, '\n');
	}
      else
	{
	  outfile = in.filename();
	  outfile.Root();
	  outfile.Append (".out");
	}
     
      ofstream log (outfile);
      log << "\n Welcome to: " << program_name << " " << program_version;
      log << "\n\n Build Date: " << __DATE__ " " << __TIME__ << "\n\n";
      MakeCSE_FEManager maker (log);

      /* read geometry and parameters */
      maker.InitializeInput (in, interactive);

      /* initialize and register to OutputBaseT */
      maker.InitializeOutput (in.filename(), program_name, program_version);

      /* make cohesive surfaces */
      maker.CreateCSE ();

      /* print output data */
      maker.WriteOutput ();
      
      cout << "\n Program Complete\n\n" << endl;
    }

  catch (int ErrorCode) 
    {
      switch (ErrorCode)
	{
	case eBadInputValue:
	  cout << "\n\n Exiting due to bad input value.\n\n";
	  break;
	case eOutOfRange:
	  cout << "\n\n Exiting due to index out of range.\n\n";
	  break;
	case eSizeMismatch:
	  cout << "\n\n Exiting due to a size mismatch.\n\n";
	  break;
	case eOutOfMemory:
	  cout << "\n\n Out of memory.\n\n";
	  break;
	default:
	  cout << "\n\n Exiting due to failure.\n\n";
	}
      cout << "      Game Over\n\n";
    }
  return 1;
}

bool IsInteractive (ifstreamT& in)
{
  StringT infile (81);
  bool file = false;
  while (file == false)
    {
      cout << "\nEnter input file name: \n" 
	   << "(\"quit\" to exit, \"nothing\" for interactive): ";
      cin.getline (infile.Pointer(), 80, '\n');
  
      if (strncmp (infile, "nothing", 7) == 0)
	return true;
      else if (strncmp (infile, "quit", 4) == 0) 
	{
	  cout << endl;
	  exit (0);
	}
      else
	{
	  in.open (infile);
	  if (in.is_open())
	    file = true;
	}
    }
  
  return false;
}

