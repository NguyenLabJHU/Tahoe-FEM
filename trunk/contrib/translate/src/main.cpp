/* $Id: main.cpp,v 1.19 2004-05-10 01:28:48 paklein Exp $ */
#include "ExceptionT.h"
#include "TranslateIOManager.h"
#include "ExtractNode.h"
#include "ExtractElement.h"
#include "ExtractQuad.h"
#include "PointPlots.h"
#include "MergeResults.h"
#include "JoinResults.h"
#include "ifstreamT.h"
#include "StringT.h"
#include "AutoArrayT.h"

using namespace Tahoe;

void ReadArgs (char* a, AutoArrayT<StringT>& list);
istream& Open (ifstreamT& tmp, const StringT& f);

int main (int c, char* a [])
{
  try 
    {
      bool write = true;
      bool echo = false;
      StringT echofile;
      AutoArrayT<StringT> filelist;
      ifstreamT tmp;

      /* is a command line arguement provided */
      if (c >= 2)
	ReadArgs (a[1], filelist);
      else
	{
	  filelist.Free();
	  StringT answer;
	  cout << "\n Echo answers to file (y/n) ? ";
	  cin >> answer;
	  if (answer[0] == 'y' || answer[0] == 'Y')
	    {
	      echo = true;
	      cout << "\n Enter file to echo parameters to: ";
	      cin >> echofile;
	    }
	}

      int numruns = 1;
      if (filelist.Length() > 0) 
	{
	  write = false;
	  numruns = filelist.Length();
	}
      for (int f=0; f < numruns; f++)
	{
	  if (write == true)
	    filelist.Append ("");
	  istream& in = Open (tmp, filelist[f]);

	  int selection;
	  if (write)
	    {
	      cout << "\n1. Datafile Translation \n";
	      cout << "2. Nodal Data Extraction to XY Data \n";
	      cout << "3. Quadrature Data Extraction to XY Data \n";
	      cout << "4. Quadrature Data Extraction for Point Plots \n";
	      cout << "5. Merge Results Files \n";
	      cout << "6. Element Data Extraction to XY Data \n";
	      cout << "7. Concatentate Results Files \n";
		  cout << "\n Select type of translation: ";
	    }
	  in >> selection;
	  cout << "\n Type of translation: " << selection << ".";

	  
	  TranslateIOManager *dataio = NULL;
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
		case 6:
		{
			cout << " Extract nodal data.\n\n";
			program = "Extract";
			version = "v1.0";
			dataio = new ExtractElement(cout, in, write);
			break;
		}
		case 7:
		{
			cout << " Concatenate results files.\n\n";
			program = "Concat";
			version = "v1.0";
			dataio = new JoinResults(cout, in, write);
			break;
		}
	    default:
	      throw ExceptionT::kGeneralFail;
	    }
	  if (echo) dataio->SetEcho (selection, echofile);
	  dataio->Translate (program, version, program);
	  delete dataio;
	}
      cout << "\n\n Progam Complete.\n\n";
    }
  catch (int ErrorCode)
    {
      cout << "\n\n Exiting due to error . . . ";
      switch (ErrorCode)
	{
	case ExceptionT::kBadInputValue:
	  cout << " Bad Input Value\n";
	  break;
	case ExceptionT::kOutOfRange:
	  cout << " Out of Range\n";
	  break;
	case ExceptionT::kSizeMismatch:
	  cout << " Size Mismatch\n";
	  break;
	case ExceptionT::kOutOfMemory:
	  cout << " Out of Memory\n";
	  break;
	case ExceptionT::kDatabaseFail:
	  cout << " Error with database\n";
	  break;
	}
      cout << "\n\n Game Over\n\n";
    }
  return 1;
}

void ReadArgs (char* a, AutoArrayT<StringT>& list)
{
  ifstreamT tmp (a);
  StringT s;
  s.GetLineFromStream (tmp);
  switch (s[0])
    {
    case '%':
      {
	list.Append (a);
	break;
      }
    case '@':
      {
	while (tmp.good())
	  {
	    s.Clear();
	    s.GetLineFromStream(tmp);
	    if (s.StringLength() > 1)
	      list.Append (s);
	  }
	break;
      }
    default:
      {
	cout << "\n You must put either % or @ at the beginning of your input file.\n\n";
	throw ExceptionT::kBadInputValue;
      }
    }
}

istream& Open (ifstreamT& tmp, const StringT& f)
{
  if (f.StringLength() >= 1)
    {
      tmp.open (f);
      cout << "\n Reading answers from: " << f << endl;

      char c;
      tmp >> c;
      if (c != '%')
	{
	  cout << "\nExpecting % at the start of the input file.\n\n";
	  throw ExceptionT::kBadInputValue;
	}
      return tmp;
    }

  return cin;
}
