// main.cpp

// created: 6 Oct 1999 by S. A. Wimmer

// program reads input file, runs MakeCSE, writes output file

#include "iosfwd.h"
#include "StringT.h"
#include "FEManager.h"
#include "MakeCSEIOManager.h"
#include "ifstream_x.h"

bool IsInteractive (ifstream_x& in);

void main (void)
{
  try
    {
      /* determine input file */
      char *program_name = "MakeCSE";
      char *program_version = "v5 (3 May 2000)";
 
      /* create input/output manager and read input file */
      cout << "\n Welcome to: " << program_name << " " << program_version;
      cout << "\n\n Build Date: " << __DATE__ " " << __TIME__ << "\n\n";

      ifstream_x in ('#');
      bool interactive = IsInteractive (in);

      StringT outfile (81);
      if (interactive)
	{
	  cout << "\n Enter output file: ";
	  cin.getline (outfile.Pointer(), 80, '\n');
	}
      else
	outfile.DefaultName (in.filename(), ".in", ".out", -1);
     
      ofstream log (outfile);
      MakeCSEIOManager iodata (log);
      iodata.ReadParameters (in, interactive, program_name, program_version);

      /* set up node and element data */
      FEManager maker (log, iodata);

      /* make cohesive surfaces */
      maker.CreateCSE ();

      /* print output data */
      maker.SetIO (iodata);
      iodata.WriteGeometry ();
      
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
}

bool IsInteractive (ifstream_x& in)
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

