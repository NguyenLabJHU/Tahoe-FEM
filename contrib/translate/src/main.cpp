/* $Id: main.cpp,v 1.10 2002-05-19 17:45:42 paklein Exp $ */

#include "TranslateIOManager.h"
#include "ExtractNode.h"
#include "ExtractQuad.h"
#include "PointPlots.h"
#include "MergeResults.h"
#include "ifstreamT.h"
#include "StringT.h"

istream& Open (int c, char* a [], ifstreamT& tmp, bool& w);

int main (int c, char* a [])
{
  try 
    {
      bool write = true;
      ifstreamT tmp;
      istream& in = Open (c, a, tmp, write);

      int selection;
      if (write)
	{
	  cout << "\n1. Datafile Translation \n";
	  cout << "2. Nodal Data Extraction to XY Data \n";
	  cout << "3. Quadrature Data Extraction to XY Data \n";
	  cout << "4. Quadrature Data Extraction for Point Plots \n";
	  cout << "5. Merge Results Files \n";
	  cout << "\n Select type of translation: ";
	}
      in >> selection;
      cout << "\n Type of translation: " << selection << ".";
      
      TranslateIOManager *dataio;
      StringT program, version;
      switch (selection)
	{
	case 1:
	  {
	    cout << " Translate data files.\n\n";
	    program = "Translate";
	    version = "v1.4";
	    dataio = new TranslateIOManager (cout, in, write);
	    break;
	  }
	case 2:
	  {
	    cout << " Extract nodal data.\n\n";
	    program = "Extract";
	    version = "v1.1";
	    dataio = new ExtractNode (cout, in, write);
	    break;
	  }
	case 3:
	  {
	    cout << " Extract quadrature data.\n\n";
	    program = "Extract";
	    version = "v1.1";
	    dataio = new ExtractQuad (cout, in, write);
	    break;
	  }
	case 4:
	  {
	    cout << " Extract quadrature data for point plots.\n\n";
	    program = "PointPlot";
	    version = "v1.0";
	    dataio = new PointPlots (cout, in, write);
	    break;
	  }
	case 5:
	  {
	    cout << " Merge data from multiple files.\n\n";
	    program = "Merge";
	    version = "v1.0";
	    dataio = new MergeResults (cout, in, write);
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

istream& Open (int c, char* a [], ifstreamT& tmp, bool& w)
{
  if (c == 2)
    {
      tmp.open (a[1]);
      cout << "\n Reading answers from: " << a[1] << endl;
      w = false;
      return tmp;
    }
  w = true;
  return cin;
}
