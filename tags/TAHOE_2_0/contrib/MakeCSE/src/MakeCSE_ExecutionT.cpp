#include "MakeCSE_ExecutionT.h"
#include "MakeCSE_FEManager.h"
#include "ExceptionT.h"
#include "ifstreamT.h"

using namespace Tahoe;

MakeCSE_ExecutionT::MakeCSE_ExecutionT (void) :
  fProgram ("MakeCSE"),
  fVersion ("v5.1 (Oct 2002)"),
  fInteractive (true)
{
  cout << "\n Welcome to: " << fProgram << " " << fVersion;
  cout << "\n\n Build Date: " << __DATE__ " " << __TIME__ << "\n\n";
}

void MakeCSE_ExecutionT::Run (const sArrayT& lineoptions)
{
  int index;
  lineoptions.HasValue ("-f", index);
  if (index > 0 && index < lineoptions.Length())
    {
      ifstreamT in ('#', lineoptions [index + 1]);
      if (!in.is_open ())
	cout << "\n Unable to open -f file: " << lineoptions [index + 1] << endl;
      else
	{
	  fInteractive = false;
	  RunBatchOrJob (in);
	}
    }
  else
    RunInteractive ();
}

/***** PRIVATE ***********/

void MakeCSE_ExecutionT::RunInteractive (void)
{
  bool done = false;
  while (!done)
    {
      StringT answer;
      cout << "\nEnter input file name: \n" 
	   << "(\"quit\" to exit, \"nothing\" for interactive): ";
      cin >> answer;
      
      if (answer == "quit")
	{
	  done = true;
	}
      else if (answer == "nothing")
	{
	  fInteractive = true;
	  ifstreamT in ('#');
	  RunJob (in);
	}
      else
	{
	  fInteractive = false;
	  ifstreamT in ('#', answer);
	  if (!in.is_open())
	    cout << "\n Unable to open: " << answer << endl;
	  else
	    RunBatchOrJob (in);
	}
    }
}

void MakeCSE_ExecutionT::RunBatchOrJob (ifstreamT& in)
{
  char filetype;
  in >> filetype;

  if (filetype == '@')
    {
      cout << "\n\n\n Reading from batch file: " << in.filename () << "\n";
      StringT nextfile;
      in >> nextfile;

      while (in.good())
	{
	  /* fix up file name */
	  //nextfile.ToNativePathName();
	  StringT path;
	  path.FilePath (in.filename());
	  nextfile.Prepend (path);

	  ifstreamT in2 ('#', nextfile);
	  if (in2.good())
	    RunBatchOrJob (in2);
	  else
	    {
	      cout << "\n Error in batchfile: " << in.filename ();
	      cout << "\n     Unable to open: " << nextfile << endl;
	      return;
	    }
	  in >> nextfile;
	}
    }
  else if (filetype == '%')
    RunJob (in);
  else
    {
      cout << "\n Unrecognized file type character: " << filetype << endl;
      return;
    }
}

void MakeCSE_ExecutionT::RunJob (ifstreamT& in)
{
      StringT outfile (81);
      if (fInteractive)
	{
	  cout << "\n Enter log file: ";
	  cin >> outfile;
	  StringT line (81);
	  cin.getline (line.Pointer(), 80, '\n'); // clear away end line char
	}
      else
	{
	  outfile = in.filename();
	  outfile.Root();
	  outfile.Append (".out");
	}
     
      cout << "\n\n Running Job: " << in.filename() << "\n\n";
      ofstream log (outfile);
      log << "\n Welcome to: " << fProgram << " " << fVersion;
      log << "\n\n Build Date: " << __DATE__ " " << __TIME__ << "\n\n";

      try
	{
	  MakeCSE_FEManager maker (log);

	  /* read geometry and parameters */
	  maker.InitializeInput (in, fInteractive);

	  /* initialize and register to OutputBaseT */
	  maker.InitializeOutput (in.filename(), fProgram, fVersion);

	  /* make cohesive surfaces */
	  maker.CreateCSE ();

	  /* print output data */
	  maker.WriteOutput ();
	}
      catch (ExceptionT::CodeT code)
	{
	  ExceptionT temp;
	  cout << "\n \"" << in.filename() << "\" exit on exception: " 
	       << code << ": " << temp.ToString(code) << "\n\n";
	  log << "\n \"" << in.filename() << "\" exit on exception: " 
	      << code << ": " << temp.ToString(code) << "\n\n";
	}
}
